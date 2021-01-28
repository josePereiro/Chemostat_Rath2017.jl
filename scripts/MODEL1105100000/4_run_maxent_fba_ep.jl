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

## ----------------------------------------------------------------------------
using Distributed

NO_WORKERS = min(length(Sys.cpu_info()) - 1, wcount)
length(workers()) < NO_WORKERS && 
    addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

## Loading everywhere
@everywhere begin

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
# GLOBAL PARAMS
@everywhere begin

    const sim_params = Dict()
    # This determine how often EP results will be cached
    sim_params[:epochlen] = 10 

    const scaling_params = Dict()
    # A smaller base could kill the process because of memory usage
    scaling_params[:scale_base] = 1000.0 
    
    const epmodel_kwargs = Dict()
    epmodel_kwargs[:alpha] = Inf

    const epconv_kwargs = Dict()
    # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:maxiter] = Int(1e4) 
    # The error threshold of convergence
    epconv_kwargs[:epsconv] = 1e-6 
    epconv_kwargs[:maxvar] = 1e50
    epconv_kwargs[:minvar] = 1e-50

    params_hash = hash((sim_params, scaling_params, epmodel_kwargs, epconv_kwargs))
end

## ----------------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "MAXENT_FBA_EP_v1"

## ----------------------------------------------------------------------------
# LOAD MODEL
UJL.tagprintln_inmw("LOADING EC MODELS")
src_file = M.FVA_PP_BASE_MODEL_FILE
fva_pp_model = UJL.load_data(src_file)
UJL.println_ifmw("size: ", size(fva_pp_model), 
    " S ChU.nzabs_range: ", ChU.nzabs_range(fva_pp_model.S), 
    "\n"
)

## ----------------------------------------------------------------------------
# CACHE MODELS
@everywhere models_cache_id = (:MODELS, sim_global_id)
UJL.save_cache(models_cache_id, fva_pp_model; headline = "MODELS CACHE SAVED")
# free 
fva_pp_model = nothing
GC.gc()

## ----------------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function prepare_model(ξ, stst)
    model = UJL.load_cache(models_cache_id; verbose = false)
    model = ChU.uncompressed_model(model)
    isnothing(model) && error("Unable to load model!!")

    # intake info
    intake_info = M.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossio's paper, (see README)
    ChSS.apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)
    return model
end

## ----------------------------------------------------------------------------
# RES IDS
# Collect all the computed results ids for bundling
const chnl = RemoteChannel() do
    Channel{Any}(10)
end
const res_ids = []
const collector = @async while true
    id = take!(chnl)
    push!(res_ids, id)
end

## ----------------------------------------------------------------------------
# SIMULATION
# Any of the loops can be parallelized by just 
# changing one of the 'map' functions
ststs = Rd.ststs
pmap(ststs) do stst
    
    ## SIMULATION PARAMS
    ξs = [Rd.val(:ξ, stst)]
    βs = [0.0; UJL.logspace(5, 7, 50)] 
    
    to_map = Iterators.product(ξs, [βs], [stst])
    map(to_map) do (ξ, βs, stst)
        
        ## HASH SEEDS
        # TODO: add PARAMS hash
        stst_hash = (stst, sim_global_id)
        sim_hash = (ξ, stst_hash)

        dat = ChSU.cached_simulation(;
            epochlen = sim_params[:epochlen], 
            verbose = true,
            sim_id = sim_hash,
            get_model = () -> prepare_model(ξ, stst),
            objider = M.BIOMASS_IDER, 
            beta_info = [(M.BIOMASS_IDER, βs)],
            costider = M.COST_IDER,
            clear_cache = false,
            use_seed = true,
            epmodel_kwargs = epmodel_kwargs,
            epconv_kwargs = epconv_kwargs
        )
        
        ## SAVING DATA
        model = prepare_model(ξ, stst)
        res_id = (:RESULT, sim_hash)
        UJL.save_cache(res_id, (stst, ξ, βs, model, dat); 
            headline = "CATCHING RESULTS\n"
        )
        
        ## PASSING ID TO MASTER
        put!(chnl, res_id)
        
        GC.gc()
        return nothing
    end # map(ξs) do ξ

    return nothing
end # map(Rd.ststs) do stst

## ----------------------------------------------------------------------------
# COLLECTING RESULTS
UJL.tagprintln_inmw("COLLECTING RESULTS ")
sleep(1) # wait for collector to get all ids
bundles = Dict()
for id in res_ids

    stst, ξ, βs, model, dat = UJL.load_cache(id; verbose = false)
    
    # Bundle
    bundle = get!(bundles, stst, ChU.ChstatBundle())

    bundle[ξ, :net] = model
    bundle[ξ, :ChLP.fba] = dat[:ChLP.fba]

    for (βi, β) in βs |> enumerate
        bundle[ξ, β, :ep] = dat[(:ep, βi)]
    end

end

## ----------------------------------------------------------------------------
# SAVING
UJL.tagprintln_inmw("SAVING RESULTS ")
UJL.save_data(M.MAXENT_FBA_EB_BUNDLES_FILE, bundles)

## ----------------------------------------------------------------------------
## CLEAR CACHE (WARNING)
if finish_clear_flag
    UJL.tagprintln_inmw("CLEARING CACHE ")
    UJL.delete_temp_caches()
end