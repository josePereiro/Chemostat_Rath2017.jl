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
NO_WORKERS = NO_CORES - 1
length(workers()) < NO_WORKERS && addprocs(NO_WORKERS; 
    exeflags = "--project")
atexit(interrupt)
println("Working in: ", workers())


# +
@everywhere begin
    
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import MAT
using SparseArrays
using StatsBase
# import JSON

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, HumanGEM, RathData, Rep_Human1, 
                            print_action, load_cached, save_cache, set_cache_dir,
                            delete_temp_caches, temp_cache_file
const HG = HumanGEM
const Rd = RathData
const RepH1 = Rep_Human1;

# A place to store temporary cached results
set_cache_dir(RepH1.MODEL_CACHE_DATA_DIR);
    
end
# -

# ---
# ### Load input

# +
@everywhere begin
    
COMP_FVA_HG_INPUT_FILE = "compareFVA_humanGEM_input.bson" # TODO: Pakage this
input_dat = wload(COMP_FVA_HG_INPUT_FILE);
orig_model = Ch.Utils.uncompress_model(input_dat["orig_model"]);
ec_model = Ch.Utils.uncompress_model(input_dat["ec_model"]);
    
end

println("\nLoaded Input at: ", relpath(COMP_FVA_HG_INPUT_FILE))
println("Orig model, size: ", size(orig_model))
println("Ec model, size: ", size(ec_model))
# -

# globals
@everywhere begin
    obj_ider = "biomass_human"
    prot_pool_exchange = "prot_pool_exchange"
    zeroth = 1e-8; # a minimum threshold to not be consider zero
end

# ---
# ### Work function

@everywhere function process_epoch(epoch, model_sym)
    
    # models are already in all workers as a global
    model = eval(model_sym)
    state = (hash(model), model_sym, epoch.start, epoch.stop)
    
    # Loading cache
    data = load_cached(state)
    !isnothing(data) && return data
    
    # Info
    print_action((model_sym, epoch.start, epoch.stop), "START EPOCH")
    
    lvals, uvals = Ch.LP.fva(model, epoch; verbose = false, on_empty_sol = (x...) -> 0.0)
    
    data = (epoch, lvals, uvals)
    save_cache(data, state)
    
    return data
end

# ---
# ### Caller function

function pmap_fva(model_sym)
    
    # models are already in all workers as a global
    model = deepcopy(eval(model_sym))
    
    T = eltype(model.S)
    M, N = size(model);
    N = 30 # Test
    epoch_len = 10 # Change here the size of each epoch
    epoch_len = clamp(epoch_len, 1, N)
    epochs = map(1:epoch_len:N) do r0
        r1 = clamp(r0 + epoch_len - 1, 0, N)
        r0:r1
    end;

    # parallel fva
    remote_res = pmap(epochs) do (epoch)
        process_epoch(epoch, model_sym)
    end;
    
    # group results
    lvals, uvals = similar(model.lb), similar(model.ub)
    for (epoch, lvals_, uvals_) in remote_res
        lvals[epoch] .= lvals_
        uvals[epoch] .= uvals_
    end
    
    model.lb .= lvals
    model.ub .= uvals
    return Ch.Utils.compress_model(model)
end

compareFVA_res = Dict()
foreach([:orig_model, :ec_model]) do model_sym
    print_action(model_sym, "STARTING FVA")
    compareFVA_res[model_sym] = pmap_fva(model_sym)
end

# Saving
# TODO: package this
COMP_FVA_HG_OUTPUT_FILE = "compareFVA_humanGEM_output.bson"
tagsave(COMP_FVA_HG_OUTPUT_FILE, Dict(DATA_KEY => compareFVA_res))
println(relpath(COMP_FVA_HG_OUTPUT_FILE), " created, size: ", filesize(COMP_FVA_HG_OUTPUT_FILE), " bytes")

delete_temp_caches()




