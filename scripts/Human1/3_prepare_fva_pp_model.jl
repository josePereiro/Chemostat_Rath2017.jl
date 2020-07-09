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
H1 = Chemostat_Rath2017.Human1
Rd = Chemostat_Rath2017.RathData
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
# -
# ### Loading everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES && addprocs(NO_CORES - 1; 
    exeflags = "--project")
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
H1 = Chemostat_Rath2017.Human1
Rd = Chemostat_Rath2017.RathData
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
    
end
# -

# This file is the primary input to the processing
if !isfile(H1.MODEL_RAW_MAT_FILE)
    error("$(H1.MODEL_RAW_MAT_FILE) not found, you must run 'make all' fisrt (see README)!!!")
end

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
    joinpath(H1.MODEL_CACHE_DATA_DIR, "$(temp_cache_file_prefix)___state_$(hash(state)).jls")

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
    
    model = deserialize(H1.BASE_MODEL_FILE);
    m, n = size(model)
    
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
model = deserialize(H1.BASE_MODEL_FILE);
m, n = size(model)
epoch_len = floor(Int, 100)
@assert epoch_len > 0

# This is parallizable
states = Iterators.product(1:epoch_len:n, epoch_len)
fva_res = pmap(process_epoch, states);

# +
# joining fva_res
lb_, ub_ = (model.lb, model.ub) .|> copy 
for (epoch, (fva_lb, fva_ub)) in fva_res
    lb_[epoch] .= fva_lb
    ub_[epoch] .= fva_ub
end

ignored = ["HMR_9136"] # put here the reactions you wants to ignore the process
ignored_idxs = [Ch.Utils.rxnindex(model, rxn) for rxn in ignored]
non_ignored = trues(n)
non_ignored[ignored_idxs] .= false

model.lb[non_ignored] = lb_[non_ignored]
model.ub[non_ignored] = ub_[non_ignored]

# deleting blocked
model = Ch.Utils.del_blocked(model; protected = ignored);
Ch.Utils.summary(model)
# -

serialize(H1.FVA_PP_BASE_MODEL_FILE, model)
println("$(relpath(H1.FVA_PP_BASE_MODEL_FILE)) created")

delete_temp_caches()


