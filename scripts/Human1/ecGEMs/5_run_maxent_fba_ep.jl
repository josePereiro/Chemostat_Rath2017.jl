## ----------------------------------------------------------------------------------------
# ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w"
        help = "number of workers to use"
        default = "1"
    "--init-clear"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--finish-clear"
        help = "clear cache at the end"   
        action = :store_true
end

if isinteractive()
    # Dev values
    wcount = 0
    init_clear_flag = false
    finish_clear_flag = false
else
    parsed_args = parse_args(set)
    wcount = parse(Int, parsed_args["w"])
    init_clear_flag = parsed_args["init-clear"]
    finish_clear_flag = parsed_args["finish-clear"]
end

# ----------------------------------------------------------------------------------------
using Distributed

NO_WORKERS = min(length(Sys.cpu_info()) - 1, wcount)
length(workers()) < NO_WORKERS && 
    addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

# ----------------------------------------------------------------------------------------
# Loading everywhere
@everywhere begin

    using DrWatson 
    quickactivate(@__DIR__, "Chemostat_Rath2017")


    using Distributed
    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    const Ch = Chemostat
    const ChU = Chemostat.Utils
    const ChSU = Chemostat.SimulationUtils
    const ChSS = Chemostat.SteadyState

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const ecG = H1.ecGEMs
    const HG = H1.HumanGEM

    import UtilsJL
    const UJL = UtilsJL

    using Plots

    UJL.set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)
    
end

## ----------------------------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    UJL.println_inmw("CLEARING CACHE ")
    UJL.delete_temp_caches()
    UJL.println_inmw("\n")
end

## ----------------------------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere FILE_ID = "5"

## ----------------------------------------------------------------------------
# GLOBALS 
@everywhere begin
    
    FIG_DIR = joinpath(ecG.MODEL_FIGURES_DATA_DIR, "$(FILE_ID)_err_progress")

    ME_BOUNDED = :ME_BOUNDED
    ME_EXPECTED = :ME_EXPECTED

end
mkpath(FIG_DIR);

## ----------------------------------------------------------------------------------------
# LOAD MODELS
UJL.println_inmw("LOADING EC MODELS")
src_file = ecG.FVA_PP_BASE_MODELS
ec_models = UJL.load_data(src_file)
model_ids = ec_models |> keys |> collect
for (model_id, model_dict) in ec_models
    model = model_dict |> ChU.compressed_copy |> ChU.MetNet
    ec_models[model_id] = model
    ChU.clampfields!(model, [:lb, :ub, :b]; abs_max = H1.MAX_BOUND, zeroth =  H1.ZEROTH)
    UJL.println_ifmw("model: ", model_id, " size: ", size(model), 
        " S ChU.nzabs_range: ", ChU.nzabs_range(model.S), "\n")
end    

## ----------------------------------------------------------------------------------------
# SCALE MODELS
UJL.println_inmw("SCALING MODELS")
scale_base = 1000.0
for (model_id, model) in ec_models
    UJL.println_ifmw("model: ", model_id, " size: ", size(model), " S ChU.nzabs_range: ", ChU.nzabs_range(model.S), "\n")
    model = ec_models[model_id] = ChU.well_scaled_model(model, scale_base; verbose = true)
    UJL.println_ifmw("ec model: ", model_id, " size: ", size(model), " S ChU.nzabs_range: ", ChU.nzabs_range(model.S), "\n")
end  

## ----------------------------------------------------------------------------------------
# CACHE MODELS
@everywhere models_cache_id = (:MODELS, FILE_ID)
UJL.save_cache(models_cache_id, ec_models; headline = "MODELS CACHE SAVED")
# free 
ec_models = nothing
GC.gc()

## ----------------------------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function load_model(model_id, stst; compressed = true)

    ξ = Rd.val(:ξ, stst)
    dat = UJL.load_cache(models_cache_id; verbose = false)
    isnothing(dat) && error("Unable to load model!!")
    model = dat[model_id] 

    # intake info
    intake_info = HG.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossio's paper, (see README)
    ChSS.apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)

    # Fix total_prot
    ChU.ub!(model, H1.PROT_POOL_EXCHANGE, 0.298) # From fba

    model = compressed ? model : ChU.uncompressed_model(model) 
    return model
end

## ----------------------------------------------------------------------------
# AUX FUNCTIONS
@everywhere function dat_file(name, ext = "jls"; kwargs...) 
    fname = UJL.mysavename(name; A = 1)
    joinpath(ecG.MODEL_PROCESSED_DATA_DIR, 
        UJL.mysavename(string(FILE_ID, "_", name), ext; kwargs...)
    )
end

@everywhere function plot_progress(sim_id, stst, method, 
        errs_tracking, max_beta_tracking
    )
    try
        errs_len = length(errs_tracking)
        curr_err = last(errs_tracking)

        UJL.tagprintln_inmw("PLOTTING", 
            "\nsim id:                      ", sim_id, 
            "\nerrs_len:                    ", errs_len, 
            "\n"
        )

        p1 = plot(log10.(errs_tracking); title = "EP error evolution",
            label = "", xlabel = "ep iteration", ylabel = "log(err)",
            lw = 3
        )
        p2 = plot(log10.(max_beta_tracking); title = "EP beta evolution",
            label = "", xlabel = "ep iteration", ylabel = "log(beta + 1e-2)",
            lw = 3
        )

        ps = Plots.Plot[p1, p2]
        layout = (2, 1)
        UJL.mysavefig(ps, "err_evolution", FIG_DIR; 
            method, layout, stst, errs_len
        )
        UJL.mysavefig(ps, "tot_err_evolution", FIG_DIR; 
            method, layout, stst
        )

    catch err; @warn("ERROR", stst, err); println() end
end

@everywhere function oniter(it, epmodel, sim_id, stst, method, 
        errs_tracking, max_beta_tracking, plot_frec
    )

    # Collect
    max_err = (it <= 1 && !isempty(errs_tracking)) ? 
        last(errs_tracking) : 
        epmodel.stat[:max_err]
    push!(errs_tracking, max_err)
    max_beta = max(maximum(epmodel.beta_vec), 1e-2)
    push!(max_beta_tracking, max_beta)
    
    errs_len = length(errs_tracking)
    try_plotting = errs_len > 1 && rem(errs_len, plot_frec) == 0.0
    try_plotting && plot_progress(sim_id, stst, method, 
        errs_tracking, max_beta_tracking
    )

    return (false, nothing)
end

## ----------------------------------------------------------------------------
INDEX_CH = RemoteChannel(() -> Channel{Any}(Inf))

## ----------------------------------------------------------------------------
# ME_EXPECTED
let
    method = ME_EXPECTED
    
    ststs = Rd.ststs[1:1] # Test
    to_iter = Iterators.product(ststs, model_ids)
    pmap(to_iter) do (stst, model_id)
        
        sim_id = (stst, method, FILE_ID)
        dfile = dat_file("me_fba_dat"; stst, method)
        put!(INDEX_CH, (stst, method, dfile))
        
        ## ----------------------------------------------------------------------------
        # CACHE
        if isfile(dfile) 
            UJL.tagprintln_inmw("CACHE FOUND (SKIPPING) ", 
               "\nsim id:        ", sim_id, 
               "\ndfile:         ", dfile,
            )
            return nothing
        end
        
        ## ----------------------------------------------------------------------------
        # SIMULATION PARAMS
        βs = [0.0; UJL.logspace(1, 14, 100)] 

        # This determine how often EP results will be cached
        epochlen = 30
        
        epconv_kwargs = Dict()
        epconv_kwargs[:maxiter] = Int(3e2)
        epconv_kwargs[:epsconv] = 1e-5
        epconv_kwargs[:damp] = 0.9
        epconv_kwargs[:maxvar] = 1e50
        epconv_kwargs[:minvar] = 1e-50

        epmodel_kwargs = Dict()
        epmodel_kwargs[:alpha] = Inf
        
        ## ----------------------------------------------------------------------------
        # TRACKING ERR
        plot_frec = 1 # Test
        errs_tracking = []
        max_beta_tracking = []
        epconv_kwargs[:oniter] = (it, epmodel) -> 
            oniter(it, epmodel, sim_id, stst, method, 
                errs_tracking, max_beta_tracking, plot_frec
            )

        ## ----------------------------------------------------------------------------
        # GET MODEL
        get_model() = load_model(model_id, stst; compressed = false)

        ## ----------------------------------------------------------------------------
        # BREAK CONDITION
        exp_μ = Rd.val(:μ, stst)
        μ_convth = 1e-2
        objidx = ChU.rxnindex(get_model(), H1.BIOMASS_IDER)
        exp_beta = 0.0

        function on_betaiter(epout, beta_vec)
            isnothing(epout) && return (false, epout)

            curr_μ = ChU.av(epout)[objidx]
            isnan(curr_μ) && error("curr_μ = NaN")
            
            err = abs(curr_μ - exp_μ)/ exp_μ
            break_ = err < μ_convth || curr_μ > exp_μ
    
            UJL.tagprintln_inmw("BEFORE EPOCH", 
                "\nsim id:                      ", sim_id, 
                "\ncurr_μ:                      ", curr_μ, 
                "\nexp_μ:                       ", exp_μ, 
                "\ngrowth err:                  ", err,
                "\nbreak:                       ", break_,
                "\n"
            )
            
            break_ && (exp_beta = maximum(beta_vec))
            return break_
        end
        
        ## ----------------------------------------------------------------------------
        # SIMULATION
        sim_dat = ChSU.cached_simulation(;
            sim_id, epochlen, get_model,
            verbose = true,
            objider = H1.BIOMASS_IDER, 
            costider = H1.PROT_POOL_EXCHANGE,
            on_betaiter,
            beta_info = [(H1.BIOMASS_IDER, βs)],
            clear_cache = false, use_seed = true,
            epmodel_kwargs, epconv_kwargs
        )

        ## ----------------------------------------------------------------------------
        # BUNDLING
        UJL.tagprintln_inmw("SAVING SIM DAT ", 
            "\nsim id:        ", sim_id, 
        )
        
        bundle = Dict()
        bundle[:exp_beta] = exp_beta
        bundle[:epout] = sim_dat[(:ep, bundle[:exp_beta])]
        bundle[:sim_dat] = sim_dat
        bundle[:fbaout] = sim_dat[:fba]
        bundle[:model] = get_model() |> ChU.compressed_model
        serialize(dfile, bundle)

        put!(INDEX_CH, (stst, method, dfile))

        return nothing
    end
end

## ----------------------------------------------------------------------------
# # ME_BOUNDED
# let
#     method = ME_BOUNDED
    
#     ststs = Rd.ststs
#     pmap(ststs) do stst
        
#         sim_id = (stst, method, FILE_ID)
#         dfile = dat_file("me_fba_dat"; stst, method)
#         put!(INDEX_CH, (stst, method, dfile))
        
#         ## ----------------------------------------------------------------------------
#         # CACHE
#         if isfile(dfile) 
#             UJL.tagprintln_inmw("CACHE FOUND (SKIPPING) ", 
#                "\nsim id:        ", sim_id, 
#             )
#             return nothing
#         end
        
#         ## ----------------------------------------------------------------------------
#         # SIMULATION PARAMS

#         # This determine how often EP results will be cached
#         epochlen = 30
        
#         epconv_kwargs = Dict()
#         epconv_kwargs[:maxiter] = Int(5e3) 
#         epconv_kwargs[:epsconv] = 1e-5
#         epconv_kwargs[:damp] = 0.9
#         epconv_kwargs[:maxvar] = 1e50
#         epconv_kwargs[:minvar] = 1e-50

#         epmodel_kwargs = Dict()
#         epmodel_kwargs[:alpha] = Inf
        
#         ## ----------------------------------------------------------------------------
#         # TRACKING ERR
#         plot_frec = 50
#         errs_tracking = []
#         max_beta_tracking = []
#         epconv_kwargs[:oniter] = (it, epmodel) -> 
#             oniter(it, epmodel, sim_id, stst, method, errs_tracking, 
#                 max_beta_tracking, plot_frec
#             )

#         ## ----------------------------------------------------------------------------
#         # GET MODEL
#         function get_model() 
#             model = load_model(:scaled_fva_pp_model, stst)
#             exp_μ = Rd.val(:μ, stst)
#             ChU.bounds!(model, M.BIOMASS_IDER, exp_μ * 0.95, exp_μ * 1.05)
#             return model
#         end

#         ## ----------------------------------------------------------------------------
#         # SIMULATION
#         sim_dat = ChSU.cached_simulation(;
#             sim_id, epochlen, get_model,
#             verbose = true,
#             objider = M.BIOMASS_IDER, 
#             costider = M.COST_IDER,
#             clear_cache = false, use_seed = true,
#             epmodel_kwargs, epconv_kwargs
#         )

#         ## ----------------------------------------------------------------------------
#         # BUNDLING
#         UJL.tagprintln_inmw("SAVING SIM DAT ", 
#             "\nsim id:        ", sim_id, 
#         )
        
#         bundle = Dict()
#         bundle[:exp_beta] = 0.0
#         bundle[:epout] = sim_dat[(:ep, bundle[:exp_beta])]
#         bundle[:fbaout] = sim_dat[:fba]
#         bundle[:model] = get_model() |> ChU.compressed_model
#         serialize(dfile, bundle)

#         put!(INDEX_CH, (stst, method, dfile))
        
#         return nothing
#     end
# end

## ----------------------------------------------------------------------------
# Collect
UJL.tagprintln_inmw("SAVING INDEX ")
INDEX = UJL.DictTree()
while isready(INDEX_CH)
    stst, method, dfile = take!(INDEX_CH)
    INDEX[:DFILE, stst, method] = relpath(dfile, ChR.PROJ_ROOT)
end
ChU.save_data(dat_file("index", "bson"), INDEX)

