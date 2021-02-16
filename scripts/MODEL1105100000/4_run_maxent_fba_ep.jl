import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w"
        help = "number of workers to use"
        default = "1"
    "--clear-at-init"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--clear-at-finish"
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
    init_clear_flag = parsed_args["clear-at-init"]
    finish_clear_flag = parsed_args["clear-at-finish"]
end

# ----------------------------------------------------------------------------
using Distributed

NO_WORKERS = min(length(Sys.cpu_info()), wcount)
length(workers()) < NO_WORKERS && 
    addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

# Loading everywhere
@time @everywhere begin

    import DrWatson: quickactivate
    quickactivate(@__DIR__, "Chemostat_Rath2017")

    using Distributed
    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSU = Ch.SimulationUtils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    
    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const M = ChR.MODEL1105100000

    import UtilsJL
    const UJL = UtilsJL

    using Plots
    
    UJL.set_cache_dir(M.MODEL_CACHE_DATA_DIR)
end

## ----------------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    UJL.tagprintln_inmw("CLEARING CACHE ")
    UJL.delete_temp_caches()
    UJL.println_inmw("\n")
end

## ----------------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere FILE_ID = "4"

## ----------------------------------------------------------------------------
# GLOBALS 
@everywhere begin
    FIG_DIR = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(FILE_ID)_err_progress")
    ME_BOUNDED = :ME_BOUNDED
end
mkpath(FIG_DIR)        

## ----------------------------------------------------------------------------
# AUX FUNCTIONS
@everywhere function dat_file(name, ext = "jls"; kwargs...) 
    fname = UJL.mysavename(name; A = 1)
    joinpath(M.MODEL_PROCESSED_DATA_DIR, 
        UJL.mysavename(string(FILE_ID, "_", name), ext; kwargs...)
    )
end

@everywhere function load_model(name, stst; compressed = false)
    MINDEX = UJL.load_data(M.MODEL_INDEX_FILE; verbose = false)
    mfile = MINDEX[stst][name]
    model = deserialize(mfile)
    compressed ? model : ChU.uncompressed_model(model)
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
# ME_BOUNDED
let
    ststs = Rd.ststs
    pmap(ststs) do stst

        method = ME_BOUNDED
        # This determine how often EP results will be cached
        epochlen = 30
        
        epconv_kwargs = Dict()
        epconv_kwargs[:maxiter] = Int(5e3) 
        epconv_kwargs[:epsconv] = 1e-5
        epconv_kwargs[:damp] = 0.9
        epconv_kwargs[:maxvar] = 1e50
        epconv_kwargs[:minvar] = 1e-50

        epmodel_kwargs = Dict()
        epmodel_kwargs[:alpha] = Inf

        sim_id = (stst, method, FILE_ID)
        
        ## ----------------------------------------------------------------------------
        # Tracking err
        plot_frec = 50
        errs_tracking = []
        max_beta_tracking = []
        epconv_kwargs[:oniter] = (it, epmodel) -> 
            oniter(it, epmodel, sim_id, stst, method, errs_tracking, max_beta_tracking, plot_frec)

        function get_model() 
            model = load_model(:scaled_fva_pp_model, stst)
            exp_μ = Rd.val(:μ, stst)
            ChU.bounds!(model, M.BIOMASS_IDER, exp_μ * 0.95, exp_μ * 1.05)
            return model
        end

        sim_dat = ChSU.cached_simulation(;
            sim_id, epochlen, get_model,
            verbose = true,
            objider = M.BIOMASS_IDER, 
            costider = M.COST_IDER,
            clear_cache = false, use_seed = true,
            epmodel_kwargs, epconv_kwargs
        )

        UJL.tagprintln_inmw("SAVING SIM DAT ", 
            "\nsim id:        ", sim_id, 
        )
        dfile = dat_file("me_fbs_dat"; stst, method)
        serialize(dfile, sim_dat)
        put!(INDEX_CH, (stst, method, dfile))
        
        return nothing
    end
end

## ----------------------------------------------------------------------------
# Collect
UJL.tagprintln_inmw("SAVING INDEX ")
INDEX = UJL.DictTree()
while isready(INDEX_CH)
    stst, method, dfile = take!(INDEX_CH)
    INDEX[:DFILE, stst, method] = relpath(dfile, ChR.PROJ_ROOT)
end
ChU.save_data(dat_file("index", "bson"), INDEX)

## ----------------------------------------------------------------------------
# # RES IDS
# # Collect all the computed results ids for bundling
# const chnl = RemoteChannel() do
#     Channel{Any}(10)
# end
# const res_ids = []
# const collector = @async while true
#     id = take!(chnl)
#     push!(res_ids, id)
# end

# ## ----------------------------------------------------------------------------
# # BETA EXPLORATION SIMULATION
# # Any of the loops can be parallelized by just 
# # changing one of the 'map' functions

# # ststs = Rd.ststs
# ststs = ["E"] # Test
# pmap(ststs) do stst
    
#     ## ----------------------------------------------------------------------------
#     ## SIMULATION PARAMS
#     βs = [0.0; UJL.logspace(1, 14, 100)] 

#     # This determine how often EP results will be cached
#     epochlen = Int(1e8)
#     epconv_kwargs = Dict()
#     # The maximum number of iteration before EP to return, even if not converged
#     epconv_kwargs[:maxiter] = Int(3e2) 
#     # The error threshold of convergence (Big at first)
#     epconv_kwargs[:epsconv] = 1e-5
#     epconv_kwargs[:damp] = 0.9
#     epconv_kwargs[:maxvar] = 1e50
#     epconv_kwargs[:minvar] = 1e-50

#     epmodel_kwargs = Dict()
#     epmodel_kwargs[:alpha] = Inf

#     # STOP CONDITION
#     comp_model = load_model(:scaled_fva_pp_model, stst; compressed = true)
#     exp_μ = Rd.val(:μ, stst)
#     μ_convth = 1e-2

#     ## ----------------------------------------------------------------------------
#     # Tracking err
#     save_frec = 50
#     errs_tracking = []
#     max_beta_tracking = []
#     FIG_DIR = joinpath(M.MODEL_FIGURES_DATA_DIR, "4_err_progress")
#     mkpath(FIG_DIR)
#     function plot_progress()
#         try
#             errs_len = length(errs_tracking)
#             curr_err = last(errs_tracking)

#             UJL.tagprintln_inmw("PLOTTING", 
#                 "\nsim id:                      ", sim_id, 
#                 "\nerrs_len:                    ", errs_len, 
#                 "\n"
#             )

#             p1 = plot(log10.(errs_tracking); title = "EP error evolution",
#                 label = "", xlabel = "ep iteration", ylabel = "log(err)",
#                 lw = 3
#             )
#             p2 = plot(log10.(max_beta_tracking); title = "EP beta evolution",
#                 label = "", xlabel = "ep iteration", ylabel = "log(beta + eps)",
#                 lw = 3
#             )

#             ps = Plots.Plot[p1, p2]
#             layout = (2, 1)
#             UJL.mysavefig(ps, "err_evolution", FIG_DIR; layout, stst, errs_len)
#             UJL.mysavefig(ps, "tot_err_evolution", FIG_DIR; layout, stst)

#         catch err; @warn("ERROR", stst, err); println() end
#     end

#     function oniter(it, epmodel)

#         # Collect
#         push!(errs_tracking, epmodel.stat[:max_err])
#         max_beta = max(maximum(epmodel.beta_vec), 1e-2)
#         push!(max_beta_tracking, max_beta)
        
#         errs_len = length(errs_tracking)
    
#         try_plotting = errs_len > 1 && rem(errs_len, save_frec) == 0.0
#         try_plotting && plot_progress()

#         return (false, nothing)
#     end
#     epconv_kwargs[:oniter] = oniter

#     ## ----------------------------------------------------------------------------
#     # BETAS
#     sim_id = (stst, FILE_ID)
#     get_model() = load_model(:scaled_fva_pp_model, stst)

#     function on_betaiter(epout)
#         isnothing(epout) && return (false, epout)
        
#         curr_μ = ChU.av(comp_model, epout, M.BIOMASS_IDER)
#         isnan(curr_μ) && error("curr_μ = NaN")
        
#         err = abs(curr_μ - exp_μ)/ exp_μ
#         break_ = err < μ_convth || curr_μ > exp_μ

#         UJL.tagprintln_inmw("BEFORE EPOCH", 
#             "\nsim id:                      ", sim_id, 
#             "\ncurr_μ:                      ", curr_μ, 
#             "\nexp_μ:                       ", exp_μ, 
#             "\ngrowth err:                  ", err,
#             "\nbreak:                         ", break_,
#             "\n"
#         )
        
#         return break_
#     end

#     sim_dat = ChSU.cached_simulation(;
#         sim_id, epochlen, 
#         verbose = true,
#         get_model, on_betaiter,
#         objider = M.BIOMASS_IDER, 
#         costider = M.COST_IDER,
#         beta_info = [(M.BIOMASS_IDER, βs)],
#         clear_cache = false, use_seed = true,
#         epmodel_kwargs, epconv_kwargs
#     )

#     serialize("test_sim_dat.jld", sim_dat)
    
#     # ## SAVING DATA
#     # res_id = (:RESULT, sim_id)
#     # UJL.save_cache(res_id, (stst, βs, sim_dat); 
#     #     headline = "CATCHING RESULTS\n"
#     # )
    
#     # ## PASSING ID TO MASTER
#     # put!(chnl, res_id)
    
#     # GC.gc()
#     return nothing
# end # map(Rd.ststs) do stst

# # ## ----------------------------------------------------------------------------
# # # COLLECTING RESULTS
# # UJL.tagprintln_inmw("COLLECTING RESULTS ")
# # sleep(1) # wait for collector to get all ids
# # bundles = Dict()
# # for id in res_ids

# #     stst, ξ, βs, model, dat = UJL.load_cache(id; verbose = false)
    
# #     # Bundle
# #     bundle = get!(bundles, stst, ChU.ChstatBundle())

# #     bundle[ξ, :net] = model
# #     bundle[ξ, :ChLP.fba] = dat[:ChLP.fba]

# #     for (βi, β) in βs |> enumerate
# #         bundle[ξ, β, :ep] = dat[(:ep, βi)]
# #     end

# # end

# # ## ----------------------------------------------------------------------------
# # # SAVING
# # 
# # UJL.save_data(M.MAXENT_FBA_EB_BUNDLES_FILE, bundles)

# # ## ----------------------------------------------------------------------------
# # ## CLEAR CACHE (WARNING)
# # if finish_clear_flag
# #     UJL.tagprintln_inmw("CLEARING CACHE ")
# #     UJL.delete_temp_caches()
# # end