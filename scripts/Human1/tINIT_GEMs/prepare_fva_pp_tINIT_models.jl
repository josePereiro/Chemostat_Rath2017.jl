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
tIG = Chemostat_Rath2017.tINIT_GEMs

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
# -
# ## Loading everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
NO_WORKERS = NO_CORES - 1
length(workers()) < NO_WORKERS && addprocs(NO_WORKERS; 
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
tIG = Chemostat_Rath2017.tINIT_GEMs
    
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
    
end

# ---
# ## Description
# ---

# This script create a new models from the ones generated in [prepare_tINIT_base_models](./prepare_tINIT_base_models.jl) (see its description for more details). It use fva to set the bounds of the model to the closest possible state. The reaction that are fixxed are deleted and modeled as a exchange/demand (except the ones protected). 

@everywhere notebook_name = "prepare_fva_pp_tINIT_models";

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
    joinpath(tIG.MODEL_CACHE_DATA_DIR, "$(temp_cache_file_prefix)___state_$(hash(state)).jls")

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
    i0, epoch_len, model_file, obj_val = state
    i0, epoch_len = (i0, epoch_len) .|> Int

    dat = deserialize(model_file);
    model = dat.metnet
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)
    m, n = size(model)

    # fix obj
    if obj_val >= 0.0
        model.lb[obj_idx] = obj_val * 0.98
        model.ub[obj_idx] = obj_val * 1.02
    end
    
    # epoch
    i1 = (i0+epoch_len > n ? n : i0+epoch_len - 1) |> Int
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

# ---
# ## FVA Preprocess
# ---
function run_fba_test(model, obj_ider)
    fbaout = Ch.LP.fba(model, obj_ider);
    println("\nFBAout summary")
    Ch.Utils.summary(model, fbaout)

    println("\nComparing with experiments")
    model = deepcopy(model)
    for stst in Rd.ststs
        println("\nStst: $stst")
        ξ = Rd.val(:ξ, stst)
        println("exp xi: $ξ")
        exp_μ = Rd.val(:μ, stst)
        println("exp growth: $exp_μ")    
        Ch.SteadyState.apply_bound!(model, ξ, HG.base_intake_info);
        fbaout_μ = Ch.LP.fba(model, obj_ider).obj_val;
        fbaout_μ < exp_μ && @warn "fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)"
        println("fba growth: $fbaout_μ")
    end
end
# -


# ## Base models

base_files = filter((s) -> startswith(s, "base_model_"), readdir(tIG.MODEL_PROCESSED_DATA_DIR))
base_files = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, base_files);
@everywhere obj_ider = "biomass_human"

for base_file in base_files
    
    # prepare epoch
    dat = deserialize(base_file);
    build_model = dat.metnet
    obj_idx = Ch.Utils.rxnindex(build_model, obj_ider)
    model_id = dat.id
    
    println("\n\n ------------------- Processing $model_id -------------------\n")
    println("\nBuild model summary")
    Ch.Utils.summary(build_model)
    fbaout = Ch.LP.fba(build_model, obj_idx);
    
    obj_val = fbaout.obj_val
    obj_val <= 0.0 && @warn "fbaout.obj_val ($(obj_val)) <= 0.0"
    m, n = size(build_model)
    epoch_len = 100 # Change here the size of each epoch
    @assert epoch_len > 0

    # This is parallizable
    println("\nParallel processing")
    states = Iterators.product(1:epoch_len:n, epoch_len, [relpath(base_file)], obj_val)
    fva_res = pmap(process_epoch, states);

    # joining fva_res
    lb_, ub_ = (build_model.lb, build_model.ub) .|> copy 
    for (epoch, (fva_lb, fva_ub)) in fva_res
        lb_[epoch] .= fva_lb
        ub_[epoch] .= fva_ub
    end

    ignored = HG.base_intake_info |> keys |> collect # put here the reactions you to protect from deleting
    ignored_idxs = [Ch.Utils.rxnindex(build_model, rxn) for rxn in ignored]
    non_ignored = trues(n)
    non_ignored[ignored_idxs] .= lb_[ignored_idxs] .!= ub_[ignored_idxs] # Only ignore the really fixed

    build_model.lb[non_ignored] = lb_[non_ignored]
    build_model.ub[non_ignored] = ub_[non_ignored]
    build_model.ib[obj_idx] = 0.0 # open obj again

    # deleting blocked
    fva_pp_model = Ch.Utils.del_blocked(build_model; protected = ignored);
    Ch.Utils.summary(fva_pp_model)
    
    println("\nFba checking")
    run_fba_test(fva_pp_model, obj_ider)

    println("\nSaving")
    file = joinpath(tIG.MODEL_PROCESSED_DATA_DIR, "fva_pp_model_$model_id.jls")
    serialize(file, (id = model_id, metnet = fva_pp_model, src = dat.src))
    println("$(relpath(file)) created")

    println("\nDeleting caches")
    delete_temp_caches()
    
    println("\n\n ------------------- Finished $model_id -------------------\n")
end


