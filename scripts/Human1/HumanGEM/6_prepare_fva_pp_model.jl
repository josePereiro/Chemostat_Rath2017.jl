# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES - 1 && addprocs(NO_CORES - 1; 
    exeflags = "--project")
atexit(interrupt)
println("Working in: ", workers())

# +
@everywhere begin
    
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

import DataFrames: DataFrame
using Dates
import Serialization: serialize, deserialize

import Chemostat
const Ch = Chemostat
import Chemostat_Rath2017: HumanGEM, RathData, print_action, temp_cache_file, 
                            save_cache, load_cached, delete_temp_caches
const HG = HumanGEM
const Rd = RathData
    
end
# -

# This file is the primary input to the processing
if !isfile(HG.MODEL_RAW_MAT_FILE)
    error("$(HG.MODEL_RAW_MAT_FILE) not found, you must run 'make all' fisrt (see README)!!!")
end

# ---
# ## Description
# ---

# This script create a new model from the one generated in [2_prepare_base_model](./2_prepare_base_model.jl) (see its description for more details). It use fva to set the bounds of the model to the closer possible state. The reaction that are fixxed are deleted and modeled as a exchange/demand (except the ones protected). 

@everywhere notebook_name = "prepare_fva_pp_model";


@everywhere begin
    obj_ider = "biomass_human"
    obj_val = Ch.LP.fba(deserialize(HG.BASE_MODEL_FILE), obj_ider).obj_val
end


# ---
# ## FVA Preprocess
# ---
# We will reduce the bounds interval of all the reactions using the results of FVA.
# If FVA for a flux returns fva_lb == fva_lb, then the flux is blocked to lb = fva_lb, ub = fva_ub
# The method allows you to set a block eps (lb = fva_lb - eps, ub = fva_ub + eps).
# We fully blocked eps = 0, for save computations in EP.

# +
@everywhere function process_epoch(state)
    
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    cached_data = load_cached(state, HG.MODEL_CACHE_DATA_DIR)
    !isnothing(cached_data) && return cached_data
    
    print_action(state, "STARTING EPOCH")
    
    # --------------------  PREPARE  --------------------  
    i0, epoch_len = state .|> Int
    
    model = deserialize(HG.BASE_MODEL_FILE);
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)
    m, n = size(model)

    # fix obj_ider
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

# +
# prepare epoch
build_model = deserialize(HG.BASE_MODEL_FILE);
obj_idx = Ch.Utils.rxnindex(build_model, obj_ider)
m, n = size(build_model)
epoch_len = floor(Int, 100)
@assert epoch_len > 0

# This is parallizable
states = Iterators.product(1:epoch_len:n, epoch_len)
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
Ch.Utils.summary(fva_pp_model)
# -

serialize(HG.FVA_PP_BASE_MODEL_FILE, fva_pp_model)
println("$(relpath(HG.FVA_PP_BASE_MODEL_FILE)) created")

delete_temp_caches(HG.MODEL_CACHE_DATA_DIR)
