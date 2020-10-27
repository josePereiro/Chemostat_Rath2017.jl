## Starting Workers
using Distributed

NO_CORES = length(Sys.cpu_info())
NO_CORES = 3
length(workers()) < NO_CORES - 1 && addprocs(NO_CORES - 1; 
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
import Chemostat.Utils: MetNet, compressed_model, uncompressed_model, clampfields!, rxnindex, del_blocked
import Chemostat.LP: fba, fva
import Chemostat_Rath2017: DATA_KEY, RathData, Human1, 
                            print_action, temp_cache_file, set_cache_dir,
                            save_cache, load_cached, delete_temp_caches
import Chemostat_Rath2017.Human1: BIOMASS_IDER, ATPM_IDER, PROT_POOL_EXCHANGE, MAX_BOUND, ZEROTH, try_fba
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
const Rd = RathData
println("Setting cache dir at: ", relpath(set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)))
    
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
        clamp_bounds!(ec_models[model_id], MAX_BOUND, ZEROTH)
    end
    myid() == 1 && println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")
end

## Finding obj_val
@everywhere begin
    println("\nComputing obj val")
    const obj_vals = Dict{String, Float64}()
    for (model_id, model) in ec_models
        obj_vals[model_id] = fba(model, BIOMASS_IDER).obj_val
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

    print_action(state, "STARTING EPOCH")

    # --------------------  PREPARE  --------------------  
    i0::Int, epoch_len::Int, model_id = state
    
    # model
    model::MetNet = deepcopy(ec_models[model_id])
    obj_idx = rxnindex(model, BIOMASS_IDER)
    m, n = size(model)

    # epoch
    i1 = i0+epoch_len > n ? n : i0+epoch_len - 1 |> Int
    epoch = i0:i1
    
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # don't trust them
    data = load_cached(state)
    if isnothing(data)
        print_action(state, "STARTING UNSAVE FVA")
        # --------------------  UNSAVE FVA  --------------------  
        data = (epoch, (model.lb[epoch], model.ub[epoch]))
        try
            data = (epoch, fva(model, epoch; verbose = false))
        catch err
            print_action(state, "ERROR DOING FVA", "Error: $err")
            err isa InterruptException && rethrow(err)
        end
        save_cache(data, state)
    end
    
    # --------------------  CHECK OBJ VAL  --------------------  
    tol = 1e-5
    backup = (model.lb[epoch], model.ub[epoch])
    model.lb[epoch] .= data[2][1]
    model.ub[epoch] .= data[2][2]
    curr_val = fba(model, BIOMASS_IDER).obj_val
    if abs(curr_val - obj_vals[model_id]) > tol 
        # --------------------  SAVE FVA  --------------------  
        print_action(state, "STARTING SAVE FVA")

        model.lb[epoch] .= backup[1]
        model.ub[epoch] .= backup[2]
        data = (epoch, (model.lb[epoch], model.ub[epoch]))
        try
            data = (epoch, fva(model, epoch; check_obj = BIOMASS_IDER, verbose = false))
        catch err
            print_action(state, "ERROR DOING FVA", "Error: $err")
            err isa InterruptException && rethrow(err)
        end
        save_cache(data, state)
    end

    # --------------------  FINISHING  --------------------  
    print_action(state, "EPOCH FINISHED")
    
    return data
    
end

## prepare epoch
fva_pp_models = Dict()
for (model_id, model) in ec_models
    build_model = deepcopy(ec_models[model_id]);
    try_fba(build_model, BIOMASS_IDER)
    obj_idx = rxnindex(build_model, BIOMASS_IDER)
    M, N = size(build_model)
    epoch_len = min(100, N)
    @assert epoch_len > 0

    # This is parallelizable
    states = Iterators.product(1:epoch_len:N, epoch_len, [model_id])
    fva_res = pmap(process_epoch, states);

    # +
    # joining fva_res
    lb_, ub_ = (build_model.lb, build_model.ub) .|> copy 
    for (epoch, (fva_lb, fva_ub)) in fva_res
        lb_[epoch] .= fva_lb
        ub_[epoch] .= fva_ub
    end

    ignored = ["HMR_9136"] # put here the reactions you wants to ignore the process
    ignored_idxs = [rxnindex(build_model, rxn) for rxn in ignored]
    non_ignored = trues(N)
    non_ignored[ignored_idxs] .= false

    build_model.lb[non_ignored] = lb_[non_ignored]
    build_model.ub[non_ignored] = ub_[non_ignored]
    build_model.lb[obj_idx] = 0.0 # open biomass again

    # deleting blocked
    fva_pp_model = del_blocked(build_model; protected = ignored);
    println("\nfva_pp_model, ", size(fva_pp_model))
    try_fba(fva_pp_model, BIOMASS_IDER)

    fva_pp_models[model_id] = compress_model(fva_pp_model)
    break;
end
# -

## Saving
# TODO: save all as compressed dicts
# file = ecG.FVA_PP_BASE_MODELS
# tagsave(file, Dict(DATA_KEY => fva_pp_models))
# println(relpath(file), " created!!!, size: ", filesize(file), " bytes")

## Clearing
delete_temp_caches()
