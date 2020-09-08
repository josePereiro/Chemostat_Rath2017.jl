using Distributed

NO_CORES = length(Sys.cpu_info())
NO_CORES = 1 #
length(workers()) < NO_CORES - 1 && addprocs(NO_CORES - 1; 
exeflags = "--project")
println("Working in: ", workers())

## Loading everywhere
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
import Chemostat.Utils: MetNet, EPModel,
                        rxnindex, metindex, compress_dict, 
                        uncompress_dict, clampfileds!, well_scaled_model,
                        ChstatBoundle, norm_abs_stoi_err, av, va, nzabs_range,
                        struct_to_dict

import Chemostat.SimulationUtils: epoch_converge_ep!, cached_simulation, set_cache_dir, 
                        tagprintln_inmw, println_inmw, tagprintln_ifmw, println_ifmw,
                        save_cache, load_cache, delete_temp_caches
import Chemostat.SteadyState: apply_bound!

# Dev
import Chemostat.LP: fba
import Chemostat.Test: toy_model
import Chemostat.Test: empty_epout
import Chemostat.MaxEntEP: maxent_ep, converge_ep!

import Chemostat_Rath2017
import Chemostat_Rath2017: DATA_KEY, Human1, RathData, load_data, save_data
import Chemostat_Rath2017.Human1: OBJ_IDER, ATPM_IDER, PROT_POOL_EXCHANGE, 
                                MAX_BOUND, ZEROTH
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
const Rd = RathData
set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)
    
end

## ------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "maxent_fba_ep_v1"

## ------------------------------------------------------------------
# PARAMS
# TODO: redesign params handling
@everywhere begin
    const params = Dict()
    params["ep_alpha"] = 1e8 # The EP inverse temperature
    params["ep_epsconv"] = 1e-5 # The error threshold of convergence
    params["ep_maxiter"] = Int(1e4) # The maximum number of iteration beforre EP to return, even if not converged
    params["ep_epoch"] = 10 # This determine how often EP results will be cached
    params["stoi_err_th"] = 3 # 
    # params["ξs"] = Rd.val(:ξ, Rd.ststs) # 
    # params["βs"] = [0.0; range(5e3, 5e4, length = 25)] 
    params["scale_base"] = 1000.0 # A smaller base could kill the process because of memory usage
end

## ------------------------------------------------------------------
# LOAD MODELS
tagprintln_inmw("LOADING EC MODELS")
src_file = ecG.FVA_PP_BASE_MODELS
ec_models = load_data(src_file)
model_ids = ec_models |> keys |> collect
for (model_id, model_dict) in ec_models
    local model = model_dict |> compress_dict |> MetNet
    ec_models[model_id] = model
    clampfileds!(model, [:lb, :ub, :b]; abs_max = MAX_BOUND, zeroth =  ZEROTH)
    println_ifmw("model: ", model_id, " size: ", size(model), " S nzabs_range: ", nzabs_range(model.S), "\n")
end    

## ------------------------------------------------------------------
# SCALE MODELS
tagprintln_inmw("SCALING MODELS")
for (model_id, model) in ec_models
    ec_models[model_id] = well_scaled_model(model, params["scale_base"]; verbose = false)
    println_ifmw("model: ", model_id, " size: ", size(model), " S nzabs_range: ", nzabs_range(model.S), "\n")
end  

## ------------------------------------------------------------------
# CACHE MODELS
@everywhere models_cache_id = (:MODELS, sim_global_id)
save_cache(models_cache_id, ec_models; headline = "MODELS CACHE SAVED")
# free 
ec_models = nothing
GC.gc()

## ------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function prepare_model(model_id, ξ, stst)
    dat = load_cache(models_cache_id; verbose = false)
    isnothing(dat) && error("Unable to load model!!")
    model = dat[model_id] 

    # intake info
    intake_info = HG.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossio's paper, (see README)
    apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)

    return model
end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
# delete_temp_caches()

## ------------------------------------------------------------------
# SIMULATION
# TODO: make it parallelizable 
# Any of the loops can be parallelized by just changing the 'map'
# for a 'pmap'.
simdata = Dict()

epmodel_kwargs = Dict()
epmodel_kwargs[:alpha] = params["ep_alpha"] # Redesign params handling

epconv_kwargs = Dict()
epconv_kwargs[:maxiter] = params["ep_maxiter"]
epconv_kwargs[:epsconv] = params["ep_epsconv"]

βs = range(1e3, 1e5, length = 25)
βs = [0.0; βs]

map(model_ids) do model_id

    tagprintln_inmw("PROCESSING MODEL ", 
        "\nid: ", model_id, 
        "\n")

    model_hash = (model_id, sim_global_id)
    simdata[model_id] = Dict()
    
    map(Rd.ststs) do stst 
        
        simdata[model_id][stst] = Dict()
        stst_hash = (stst, model_hash)
        ξs = [Rd.val(:ξ, stst)]
        
        map(ξs) do ξ
            sim_hash = (stst, model_hash)
            simdata[model_id][stst][ξ] = sim_hash

            dat = cached_simulation(;
                epochlen = 5, # TODO: handle better
                verbose = true,
                sim_id = sim_hash,
                get_model = function()
                    return prepare_model(model_id, ξ, stst);
                end,
                objider = OBJ_IDER,
                costider = PROT_POOL_EXCHANGE,
                beta_info = [(OBJ_IDER, βs)], # TODO: handle betas better
                clear_cache = false,
                use_seed = true,
                epmodel_kwargs = epmodel_kwargs,
                epconv_kwargs = epconv_kwargs
            )
            
            save_cache(sim_hash, dat, 
                headline = "SIMULATION DATA SAVED")
            
            GC.gc()
            return nothing
        end # map(ξs) do ξ

        return nothing
    end # map(Rd.ststs) do stst

    return nothing
end # map(model_ids) do model_id

## COLLECTING RESULTS
tagprintln_inmw("COLLECTING RESULTS ")
for (model_id, model_dat) in simdata
    for (stst, stst_dat) in model_dat
        for (ξ, cache_hash) in stst_dat
            simdata[model_id][stst][ξ] = load_cache(cache_hash; verbose = false)
        end 
    end
end

## SAVING
# TODO: package this
tagprintln_inmw("SAVING RESULTS ")
save_data(ecG.MAXENT_FBA_EB_BOUNDLES_FILE, simdata)

## CLEAR CACHE (WARNING)
tagprintln_inmw("CLEARING CACHE ")
# delete_temp_caches()