# -*- coding: utf-8 -*-
# ### Precompaling in master worker first

# +
import DataFrames: DataFrame
import MAT
import CSV
using Distributed
using Dates
import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
Rd = Chemostat_Rath2017.RathData
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
# -
# ### Loading everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES - 1 && addprocs(NO_CORES - 1; 
    exeflags = "--project")
atexit(interrupt)
println("Working in: ", workers())

# +
@everywhere begin
    
import DataFrames: DataFrame
import MAT
import CSV
using Dates
using Distributed
import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
Rd = Chemostat_Rath2017.RathData
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
    
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

@everywhere notebook_name = "2_prepare_fva_pp_model";

# ---
# ## Print functions
# ---

@everywhere function print_action(wid, pid, state, head, bodyls...)
    Core.println("Worker $wid ($pid) $head at $(Time(now())) ----------------------------")
    Core.println("\tState: ", state)
    for body in bodyls
        Core.println("\t$body")
    end
    Core.println()
    flush(stdout);
end
@everywhere function print_action(state, head, bodyls...)
    remotecall_wait(print_action, 1, myid(), getpid(), state, head, bodyls...)
end

# ---
# ## temp caching
# ---

@everywhere temp_cache_file_prefix = "$(notebook_name)___temp_cache"
@everywhere temp_cache_file(state...) = 
    joinpath(HG.MODEL_CACHE_DATA_DIR, "$(temp_cache_file_prefix)___state_$(hash(state)).jls")

@everywhere function load_cached(state)
    
    tcache_file = temp_cache_file(state) |> relpath
    data = nothing
    if isfile(tcache_file)
        try
            data = deserialize(tcache_file)
        catch err
            print_action(state, "ERROR LOADING CACHE", 
                "cache_file: $tcache_file", 
                "err:        $(err)")
        end
        print_action(state, "CACHE LOADED", "cache_file: $tcache_file")
    end
    return data
end

@everywhere function save_cache(data, state)
    tcache_file = temp_cache_file(state) |> relpath
    try
        serialize(tcache_file, data)
    catch err
         print_action(state, "ERROR SAVING CACHE", 
                "cache_file: $tcache_file", 
                "err:        $(err)")
    end
    print_action(state, "CACHE SAVED", "cache_file: $tcache_file")
end

@everywhere function delete_temp_caches()
    cache_dir = dirname(temp_cache_file())
    tcaches = filter(file -> occursin(temp_cache_file_prefix, file), readdir(cache_dir))
    for tc in tcaches
        tc = joinpath(cache_dir, tc)
        rm(tc, force = true)
        println(relpath(tc), " deleted!!!")
    end
end

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
    cached_data = load_cached(state)
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
    save_cache(data, state)
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

delete_temp_caches()
