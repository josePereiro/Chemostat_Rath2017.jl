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
# GLOBAL PARAMS
@everywhere begin

    const sim_params = Dict()
    sim_params[:epochlen] = 10 # This determine how often EP results will be cached

    const scaling_params = Dict()
    scaling_params[:scale_base] = 1000.0 # A smaller base could kill the process because of memory usage
    
    const epmodel_kwargs = Dict()
    epmodel_kwargs[:alpha] = Inf

    const epconv_kwargs = Dict()
    epconv_kwargs[:maxiter] = Int(1e4) # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:epsconv] = 1e-5 # The error threshold of convergence
    epconv_kwargs[:maxvar] = 1e35
    epconv_kwargs[:minvar] = 1e-35

    params_hash = hash((sim_params, scaling_params, epmodel_kwargs, epconv_kwargs))

end

## ----------------------------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "MAXENT_FBA_EP_v2"

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
# UJL.println_inmw("SCALING MODELS")
# for (model_id, model) in ec_models
#     ec_models[model_id] = well_scaled_model(model, scaling_params[:scale_base]; verbose = false)
#     UJL.println_ifmw("model: ", model_id, " size: ", size(model), " S ChU.nzabs_range: ", ChU.nzabs_range(model.S), "\n")
# end  

## ----------------------------------------------------------------------------------------
# CACHE MODELS
@everywhere models_cache_id = (:MODELS, sim_global_id)
UJL.save_cache(models_cache_id, ec_models; headline = "MODELS CACHE SAVED")
# free 
ec_models = nothing
GC.gc()

## ----------------------------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function prepare_model(model_id, ξ, stst)
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

    return model |> ChU.uncompressed_model
end

## ----------------------------------------------------------------------------------------
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

## ----------------------------------------------------------------------------------------
# SIMULATION
# Any of the loops can be parallelized by just changing one of the 'map' functions
to_map = Iterators.product(model_ids)
map(model_ids) do (model_id)

    UJL.println_inmw("PROCESSING MODEL ", 
        "\nid: ", model_id, 
        "\n")
    
    ststs = Rd.ststs[1:1] # Test
    to_map = Iterators.product(ststs, [model_id])
    pmap(to_map) do (stst, model_id)
        
        ## SIMULATION PARAMS
        ξs = [Rd.val(:ξ, stst)]
        βs = [0.0; range(1e3, 1e5, length = 25)] 
        
        to_map = Iterators.product(ξs, [βs], [stst], [model_id])
        map(to_map) do (ξ, βs, stst, model_id)
            
            ## HASH SEEDS
            # TODO: add PARAMS hash
            model_hash = (model_id, sim_global_id)
            stst_hash = (stst, model_hash)
            sim_hash = (ξ, stst_hash)

            dat = ChSU.cached_simulation(;
                epochlen = sim_params[:epochlen], 
                # epochlen = 100, # Test
                verbose = true,
                sim_id = sim_hash,
                get_model = function()
                    return prepare_model(model_id, ξ, stst);
                end,
                objider = H1.BIOMASS_IDER, 
                beta_info = [(H1.BIOMASS_IDER, βs)],
                costider = H1.PROT_POOL_EXCHANGE,
                clear_cache = false,
                use_seed = true,
                epmodel_kwargs = epmodel_kwargs,
                epconv_kwargs = epconv_kwargs
            )
            
            ## SAVING DATA
            model = prepare_model(model_id, ξ, stst)
            res_id = (:RESULT, sim_hash)
            UJL.save_cache(res_id, (model_id, stst, ξ, βs, model, dat); 
                headline = "CATCHING RESULTS\n"
            )
            
            ## PASSING ID TO MASTER
            put!(chnl, res_id)
            
            GC.gc()
            return nothing
        end # map(ξs) do ξ

        return nothing
    end # map(Rd.ststs) do stst

    return nothing
end # map(model_ids) do model_id

## ----------------------------------------------------------------------------------------
# ## COLLECTING RESULTS
# UJL.println_inmw("COLLECTING RESULTS ")
# sleep(1) # wait for collector to get all ids
# bundles = Dict()
# for id in res_ids

#     model_id, stst, ξ, βs, model, dat = UJL.load_cache(id; verbose = false)
    
#     # Bundle
#     model_dict = get!(bundles, model_id, Dict())
#     bundle = get!(model_dict, stst, ChstatBundle())

#     bundle[ξ, :net] = model
#     bundle[ξ, :fba] = dat[:fba]

#     for (βi, β) in βs |> enumerate
#         bundle[ξ, β, :ep] = dat[(:ep, βi)]
#     end

# end

# ## ----------------------------------------------------------------------------------------
# # SAVING
# UJL.println_inmw("SAVING RESULTS ")
# save_data(ecG.MAXENT_FBA_EB_BOUNDLES_FILE, bundles)

# ## ----------------------------------------------------------------------------------------
# # CLEAR CACHE (WARNING)
# if finish_clear_flag
#     UJL.println_inmw("CLEARING CACHE ")
#     UJL.delete_temp_caches()
# end