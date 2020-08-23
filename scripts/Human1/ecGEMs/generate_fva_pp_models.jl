## Starting Workers
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES - 2 && addprocs(NO_CORES - 2; 
    exeflags = "--project")
atexit(interrupt)
println("Working in: ", workers())

## Importing everywhere
@everywhere begin
    
import DrWatson: quickactivate, wload, tagsave
quickactivate(@__DIR__, "Chemostat_Rath2017")

using Dates
using SparseArrays
import Serialization: serialize, deserialize

import Chemostat
const Ch = Chemostat
import Chemostat.Utils: MetNet, uncompress_model
import Chemostat_Rath2017: DATA_KEY, RathData, Human1, 
                            print_action, temp_cache_file, set_cache_dir,
                            save_cache, load_cached, delete_temp_caches
import Chemostat_Rath2017.Human1: OBJ_IDER, ATPM_IDER, PROT_POOL_EXCHANGE
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
const Rd = RathData
set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)
    
end

## Description
# This script create a new models from the one generated in 
# [generate_base_brain_ecModels](./generate_base_brain_ecModels.jl) (see its description for more details). 
# It use fva to set the bounds of the model to the closer possible state. 
# The reaction that are fixxed are deleted and modeled as a exchange/demand (except the protected ones). 

## Loading base models
@everywhere begin
    println("\nLoading ec base models")
    src_file = ecG.EC_BRAIN_BASE_MODELS_FILE
    const ec_models = wload(src_file)[DATA_KEY]
    for (model_id, model) in ec_models
        ec_models[model_id] = uncompress_model(model)
    end
    myid() == 1 && println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")
end

## Finding obj_val
@everywhere begin
    println("\nComputing obj val")
    const obj_vals = Dict{String, Float64}()
    for (model_id, model) in ec_models
        obj_vals[model_id] = Ch.LP.fba(model, OBJ_IDER).obj_val
        myid() == 1 && println("model ", model_id, " obj_val: ", obj_vals[model_id])
    end
end


## FVA Preprocess
# We will reduce the bounds interval of all the reactions using the results of FVA.
# If FVA for a flux returns fva_lb == fva_lb, then the flux is blocked to lb = fva_lb, ub = fva_ub
# The method allows you to set a block eps (lb = fva_lb - eps, ub = fva_ub + eps).
# We fully blocked eps = 0, for save computations in EP.

##
@everywhere function process_epoch(state)
    
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # don't trust them
    cached_data = load_cached(state, HG.MODEL_CACHE_DATA_DIR)
    !isnothing(cached_data) && return cached_data
    
    print_action(state, "STARTING EPOCH")
    
    # --------------------  PREPARE  --------------------  
    i0::Int, epoch_len::Int, model_id = state
    
    model::MetNet = deepcopy(ec_models[model_id])
    obj_idx = Ch.Utils.rxnindex(model, OBJ_IDER)
    m, n = size(model)

    # fix obj_ider
    obj_val = obj_vals[model_id]
    model.lb[obj_idx] = obj_val * 0.98
    model.ub[obj_idx] = obj_val * 1.02

    # epoch
    i1 = i0+epoch_len > n ? n : i0+epoch_len - 1 |> Int
    epoch = i0:i1

    # --------------------  FVA  --------------------  
    data = (epoch, (model.lb[epoch], model.ub[epoch]))
    try
        data = (epoch, Ch.LP.fva(model, epoch; verbose = false))
    catch err
        print_action(state, "ERROR DOING FVA", "Error: $err")
        err isa InterruptException && rethrow(err)
    end
    
    # --------------------  FINISHING  --------------------  
    save_cache(data, state, HG.MODEL_CACHE_DATA_DIR)
    print_action(state, "EPOCH FINISHED")
    
    return data
    
end

## prepare epoch
fva_pp_models = Dict()
for (model_id, model) in ec_models
    build_model = deepcopy(ec_models[model_id]);
    obj_idx = Ch.Utils.rxnindex(build_model, OBJ_IDER)
    m, n = size(build_model)
    epoch_len = floor(Int, 100)
    @assert epoch_len > 0

    # This is parallelizable
    states = Iterators.product(1:epoch_len:n, epoch_len, [model_id])
    fva_res = pmap(process_epoch, states);

    # +
    # joining fva_res
    lb_, ub_ = (build_model.lb, build_model.ub) .|> copy 
    for (epoch, (fva_lb, fva_ub)) in fva_res
        lb_[epoch] .= fva_lb
        ub_[epoch] .= fva_ub
    end

    ignored = ["HMR_9136"] # put here the reactions you wants to ignore the process
    ignored_idxs = [Ch.Utils.rxnindex(build_model, rxn) for rxn in ignored]
    non_ignored = trues(n)
    non_ignored[ignored_idxs] .= false

    build_model.lb[non_ignored] = lb_[non_ignored]
    build_model.ub[non_ignored] = ub_[non_ignored]
    build_model.lb[obj_idx] = 0.0 # open biomass again

    # deleting blocked
    fva_pp_model = Ch.Utils.del_blocked(build_model; protected = ignored);
    println("\nfva_pp_model, ", size(fva_pp_model))
    Human1.try_fba(fva_pp_model, OBJ_IDER)

    fva_pp_models[model_id] = fva_pp_model
    break;
end
# -

file = ecG.FVA_PP_BASE_MODELS
tagsave(file, Dict(DATA_KEY => fva_pp_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")

# delete_temp_caches(HG.MODEL_CACHE_DATA_DIR)
