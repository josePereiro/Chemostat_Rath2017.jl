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

using DataFrames
using Serialization
using Dates
import StatsBase: mean

# ### Precompaling in master worker first

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
HG = Chemostat_Rath2017.HumanGEM
tIG = Chemostat_Rath2017.tINIT_GEMs

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# ### Loading everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
# length(workers()) < NO_CORES - 2 && addprocs(NO_CORES - 2; 
#     exeflags = "--project")
println("Working in: ", workers())

# +
@everywhere begin

using DataFrames
using Distributed
using Serialization
using Dates

import Chemostat
Ch = Chemostat

# TODO include installation help
import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
HG = Chemostat_Rath2017.HumanGEM
tIG = Chemostat_Rath2017.tINIT_GEMs
    
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env()
    
end
# -

# ---
# ## Description
# ---

# This script use the model produced in [3_prepare_fva_pp_model](./3_prepare_fva_pp_model.jl) (see 3_prepare_fva_pp_model description for more details). It use MaxEnt EP algorithm as reported in [Cossio et al](https://doi.org/10.1371/journal.pcbi.1006823). 

@everywhere notebook_name = "fva_pp_tINIT_models_maxent_ep";

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

@everywhere function print_progress(wid, pid, 
            stst, ξi, ξs, ξ,  βi, βs, β, 
            exp_av, fba_av, 
            ep_av, ep_alpha, ep_epsconv, ep_maxiter, ep_status, ep_iter,
            min_err, mean_err, max_err,
            elapsed)
    Core.println("Worker: $wid ($pid) stst: $stst progress at $(Time(now())) ----------------------------")
    Core.println("\txi: [$ξi/ $(length(ξs))] beta: [$βi/ $(length(βs))]")
    Core.println("\t  ----------------- --------------")
    Core.println("\tmodel xi:           $ξ")
    Core.println("\tmodel beta:         $β")
    Core.println("\tmodel ep_alpha:     $ep_alpha")
    Core.println("\tmodel ep_epsconv:   $ep_epsconv")
    Core.println("\tmodel ep_status:    $ep_status")
    Core.println("\tmodel ep_iter:      $ep_iter")
    Core.println("\tmodel ep_maxiter:   $ep_maxiter")
    Core.println("\tmodel fba obj:      $fba_av")
    Core.println("\tmodel ep obj:       $ep_av")
    Core.println("\tmodel ep stoi err:  $((min_err, mean_err, max_err))")
    Core.println("\t  ----------------- --------------")
    Core.println("\texp   xi:           $(Rd.val(:ξ, stst))")
    Core.println("\texp   cGLC:         $(Rd.val(:cGLC, stst))")
    Core.println("\texp   cGLN:         $(Rd.val(:cGLN, stst))")
    Core.println("\texp   cGAL:         $(Rd.val(:cGAL, stst))")
    Core.println("\texp   obj:          $exp_av")
    Core.println("\t  ----------------- --------------")
    Core.println("\txi elapsed time(s): $elapsed")
    Core.println()
    flush(stdout);
end
@everywhere print_progress(stst, ξi, ξs, ξ,  βi, βs, β, exp_av, fba_av, ep_av, ep_alpha, 
        ep_epsconv, ep_maxiter, ep_status, ep_iter, min_err, mean_err, max_err, elapsed) = 
    remotecall_wait(print_progress, 1, myid(), getpid(), 
        stst, ξi, ξs, ξ,  βi, βs, β, exp_av, fba_av, ep_av, 
        ep_alpha, ep_epsconv, ep_maxiter, ep_status, ep_iter, min_err, mean_err, max_err, elapsed)

# ---
# ## temp caching
# ---

@everywhere temp_cache_file_prefix = "$(notebook_name)___temp_cache"
@everywhere temp_cache_file(state) = 
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

@everywhere function save_cache(state, data)
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
    cache_dir = dirname(temp_cache_file("bla"))
    tcaches = filter(file -> occursin(temp_cache_file_prefix, file), readdir(cache_dir))
    for tc in tcaches
        tc = joinpath(cache_dir, tc)
        rm(tc, force = true)
        println(relpath(tc), " deleted!!!")
    end
end

# ---
# ## Testing
# ---

@everywhere begin
    testing = false # Change this to true for testing
end
println("Testing: ", testing)

# ---
# ## Params
# ---

# +
@everywhere begin
    
    params = Dict()
    
    # model iders
    params["obj_ider"] = "biomass_human"

    # ep
    params["ep_alpha"] = 1e7 # The EP inverse temperature
    params["ep_epsconv"] = 1e-5 # The error threshold of convergence
    params["ep_maxiter"] = 1e4 # The maximum number of iteration beforre EP to return, even if not converged
    params["ep_epoch"] = 10 # This determine how often EP results will be cached

end
# -

# ---
# ## work functions
# ---

# +
@everywhere function process_exp(stst, model_file, upfrec = 10)
    
    # --------------------  CHEMOSTAT PARAMETER XI --------------------  
    # Change here how many xi to model, you should always include the experimental xis
    ξs = [Rd.val("ξ", stst)]
    ξs = testing ? [Rd.val("ξ", stst)] : ξs
    
    
    # --------------------  MAXENT PARAMETER BETA --------------------  
    # Change here how many betas to model
    # The beta range is set up by trial and error
    βs = range(5e3, 5e4, length = 25)
    βs = [0.0; βs]
    βs = testing ? collect(1:5) : βs

    # Current state
    state = (stst, model_file, hash((ξs, βs)))
    
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    cached_data = load_cached(state)
    !isnothing(cached_data) && return cached_data
    
    
    # Print hello in worker 1
    print_action(state, "STARTING EXPERIMENT")
    
    
    # --------------------  PROCESSING EACH XI --------------------  
    # Here parallel processing can optionally be used by calling 'pmap' or
    # 'map' instead, depending in the configuration of the parameters and the number
    # of experiments
    ixs_data = map((ξi) -> process_xi(stst, model_file, ξi, ξs, βs, upfrec), eachindex(ξs))
    
    # --------------------  BOUNDLING --------------------  
    boundle = Ch.Utils.ChstatBoundle()
    foreach((ξi) -> boundle_xi_data!(boundle, ξs[ξi], βs, ixs_data[ξi]), eachindex(ξs))

    # --------------------  FINISHING --------------------   
    data = (stst, boundle)
    save_cache(state, data)
    
    # Printing finishing in 1
    print_action(state, "FINISHING EXPERIMENT")
    
    return data
    
end

# +
@everywhere function process_xi(stst, model_file, ξi, ξs, βs, upfrec)
    
    # Current state
    ξ = ξs[ξi]
    xi_state = (stst, ξ, hash(βs))
    t0 = time() # to track xi processing duration
        
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    cached_data = load_cached(xi_state)
    !isnothing(cached_data) && return cached_data

    # Print hello in worker 1
    print_action(xi_state, "STARTING XI")
    
    # --------------------  EP PARAMS  --------------------  
    ep_alpha = params["ep_alpha"]
    ep_epsconv = params["ep_epsconv"]
    ep_maxiter = floor(Int, params["ep_maxiter"])
    ep_epoch = params["ep_epoch"]


    # --------------------  PREPARING MODEL  --------------------  
    dat = deserialize(model_file);
    model = dat.metnet
    m, n = size(model)
    obj_ider = params["obj_ider"]
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)
    
    # intake info
    intake_info = Dict()
    for (intake, info) in HG.stst_base_intake_info(stst)
        if intake in model.rxns 
            intake_info[intake] = info
        else
            delete!(model.intake_info, intake)
        end
    end

    # Chemostat steady state constraint, see Cossios paper, (see README)
    Ch.SteadyState.apply_bound!(model, ξ, intake_info)

    # --------------------  FBA  --------------------  
    fbaout = Ch.LP.fba(model, obj_ider)

    # storing
    xi_data = Dict() # local store
    xi_data[(ξ, :net)] = model
    xi_data[(ξ, :fba)] = fbaout
    
    
    # --------------------  BETA LOOP  --------------------  
    # THIS CAN TAKE A WHILE!!!
    βv = zeros(n)
    show_t = time()
    beta_seed = nothing
    for (βi, β) in βs |> enumerate
        
        # --------------------  TEMP CACHE  --------------------  
        # I do not check any cache consistency, so delete the temporal caches if you
        # dont trust them
        # currrent state
        beta_state = (stst, ξ, β)
        cached = load_cached(beta_state)
        epout = isnothing(cached) ? beta_seed : cached
        
        # --------------------  TRY MAXENT-EP  --------------------  
        try
            
            # --------------------  MAXENT-EP WHILE LOOP  --------------------  
            # I break down ep in several epoch to be able to cache partial results
            βv[obj_idx] = β
            curr_iter = 0
            while (ep_maxiter > curr_iter) && # Till maxiter 
                (isnothing(epout) || epout.status != :converged) # Or converged
                
                epout = testing ? Ch.Test.empty_epout(model) :
                    Ch.MaxEntEP.maxent_ep(model, 
                            alpha = ep_alpha, 
                            beta_vec = βv, 
                            epsconv = ep_epsconv, 
                            maxiter = ep_epoch,
                            solution = epout, 
                            verbose = false)
                
                curr_iter += epout.iter

                # Show progress # TODO epout.iter > 10, make this setable
                if βi == 1 || βi == length(βs) || βi % upfrec == 0 || 
                        epout.iter > 10 || epout.status == :converged ||
                        (time() - show_t) > 30 # update if more than 30s
                    show_t = time()
                    
                    exp_av = Rd.val(:μ, stst)
                    fba_av = Ch.Utils.av(model, fbaout, obj_idx)
                    ep_av = Ch.Utils.av(model, epout, obj_idx)
                    ep_status = epout.status
                    ep_iter = curr_iter

                    abs_norm_err = Ch.Utils.norm_abs_stoi_err(model, epout) .|> abs
                    max_err = abs_norm_err |> maximum
                    min_err = abs_norm_err |> minimum
                    mean_err = abs_norm_err |> mean

                    elapsed = time() - t0

                    print_progress(stst, ξi, ξs, ξ, βi, βs, β, exp_av, fba_av, ep_av, 
                        ep_alpha, ep_epsconv, ep_maxiter, ep_status, ep_iter, 
                        min_err, mean_err, max_err,
                        elapsed)
                end
                
                # caching
                save_cache(beta_state, epout)
                
            end # EP While loop

            # save for the next beta
            beta_seed = epout

            # storing
            xi_data[(ξ, β, :ep)] = epout

        catch err
            # Print error in 1
            print_action(beta_state, "ERROR DURING EP", "err: $(err)")

            # storing error
            xi_data[(ξ, β, :ep)] = err
            
            err isa InterruptException && rethrow(err)
        end # try
        
    end # β loop
    
    # --------------------  FINISHING --------------------   
    save_cache(xi_state, xi_data)
    
    # Printing finishing in 1
    print_action(xi_state, "FINISHING XI")
    
    return xi_data
    
end # process_xi
# -

@everywhere function boundle_xi_data!(boundle, ξ, βs, xi_data)

    boundle[ξ, :net] = get(xi_data, (ξ, :net), nothing)
    boundle[ξ, :fba] = get(xi_data, (ξ, :fba), nothing)
    
    for β in βs
        boundle[ξ, β, :ep] = get(xi_data, (ξ, β, :ep), nothing)
    end
end 

# ### Processing

fva_pp_files = filter((s) -> startswith(s, "fva_pp_model_"), readdir(tIG.MODEL_PROCESSED_DATA_DIR))
fva_pp_files = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, fva_pp_files);

for model_file in fva_pp_files

    model_file = relpath(model_file)
    dat = deserialize(model_file)
    model_id = dat.id
    boundle_file = joinpath(tIG.MODEL_PROCESSED_DATA_DIR, "$(notebook_name)___$(model_id)___boundles.jls")
    if isfile(boundle_file)
        println("\n skipping $model_id, $(basename(model_file)) exist!!!\n")
        continue
    end

    println("\n\n ------------------ $model_id ------------------\n\n")
    Ch.Utils.summary(dat.metnet)
    println()

    # this can take a while!!!
    # A, B, C steady states have the same initial conditions
    # ststs_ = [stst for stst in Rd.ststs if stst != "B" && stst != "C"]
    ststs_ = testing ? Rd.ststs[1:1] : Rd.ststs
    println("Ststs: ", ststs_)

    remote_results = map((stst) -> process_exp(stst, model_file), ststs_);

    # ### Saving

    println("\nSaving")
    to_save = (params = params, id = model_id, res = remote_results)
    !testing && serialize(boundle_file, to_save)
    println(relpath(boundle_file), " created!!!")

    # ### Delete Temp Caches

    # Do not forget to run this if you change any parameter
    println("\nDeleting caches")
    delete_temp_caches()
    println()
    flush(stdout)
    
end