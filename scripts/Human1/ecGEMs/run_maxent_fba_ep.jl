
using Distributed

NO_CORES = length(Sys.cpu_info())
NO_CORES = 3 #
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
import Chemostat.Utils: MetNet, rxnindex, metindex, compress_model, 
                        uncompress_model, clamp_bounds!,
                        ChstatBoundle, norm_abs_stoi_err, av, va
import Chemostat.LP: fba
import Chemostat.SteadyState: apply_bound!
import Chemostat.Test: empty_epout
import Chemostat.MaxEntEP: maxent_ep

import Chemostat_Rath2017
import Chemostat_Rath2017: DATA_KEY, Human1, print_action, string_err,
                            load_cached, save_cache, set_cache_dir,
                            delete_temp_caches, temp_cache_file, RathData,
                            println_if1
import Chemostat_Rath2017.Human1: OBJ_IDER, ATPM_IDER, PROT_POOL_EXCHANGE, 
                                MAX_BOUND, ZEROTH
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
const Rd = RathData
set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)
    
end

## Description
# This script use the model produced in [generate_fva_pp_model](./3_prepare_fva_pp_model.jl) 
# (see 3_prepare_fva_pp_model description for more details). It use MaxEnt EP algorithm as reported 
# in [Cossio et al](https://doi.org/10.1371/journal.pcbi.1006823). 

## Print functions
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

## Testing
@everywhere begin
    testing = false # make true for testing
end
println("Testing: ", testing)

## Params
@everywhere begin
    const params = Dict()

    # ep
    params["ep_alpha"] = 1e7 # The EP inverse temperature
    params["ep_epsconv"] = 1e-5 # The error threshold of convergence
    params["ep_maxiter"] = 1e4 # The maximum number of iteration beforre EP to return, even if not converged
    params["ep_epoch"] = 10 # This determine how often EP results will be cached
    params["stoi_err_th"] = 3 # 
    params["ξs"] = Rd.val(:ξ, Rd.ststs) # 
    params["βs"] = [0.0; range(5e3, 5e4, length = 25)] # 

end

## Loading models
println("\n\n ------------------ LOADING EC MODELS ------------------\n\n")
@everywhere begin
    src_file = ecG.FVA_PP_BASE_MODELS
    const ec_models = wload(src_file)[DATA_KEY]
    println_if1(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")
    for (model_id, model) in ec_models
        ec_models[model_id] = uncompress_model(model)
        clamp_bounds!(ec_models[model_id], MAX_BOUND, ZEROTH)
        println_if1("model: ", model_id, " size: ", size(model))
    end
end


## work functions
# Returns a boundle with all the data
@everywhere function process_exp(state)

    # State
    model_id, stst = state

    # --------------------  CHEMOSTAT PARAMETER XI --------------------  
    # Change here how many xi to model, you should always include the experimental xis
    ξs = [Rd.val(:ξ, stst)]
    ξs = testing ? [Rd.val("ξ", stst)] : ξs # params["ξs"]
    
    
    # --------------------  MAXENT PARAMETER BETA --------------------  
    # Change here how many betas to model
    # The beta range is set up by trial and error
    βs::Vector{Float64} = testing ? collect(1:5) : params["βs"]
    
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # don't trust them
    cache_state = (stst, model_id, hash(params))
    cached_data = load_cached(cache_state)
    !isnothing(cached_data) && return cached_data
    
    
    # Print hello in worker 1
    print_action(state, "STARTING EXPERIMENT")
    
    
    # --------------------  PROCESSING EACH XI --------------------  
    # Here parallel processing can optionally be used by calling 'pmap' or
    # 'map' instead, depending in the configuration of the parameters and the number
    # of experiments
    ixs_data = map(eachindex(ξs)) do ξi
        xi_state = (stst, model_id, ξi, ξs[ξi])
        process_xi(xi_state)
    end
    
    # --------------------  BOUNDLING --------------------  
    boundle = ChstatBoundle()
    foreach((ξi) -> boundle_xi_data!(boundle, ξs[ξi], βs, ixs_data[ξi]), eachindex(ξs))

    # --------------------  FINISHING --------------------   
    data = boundle
    save_cache(data, cache_state)
    
    # Printing finishing in 1
    print_action(state, "FINISHING EXPERIMENT")
    
    return data
    
end

# +
@everywhere function process_xi(xi_state; upfrec = 10)
    
    # Current state
    stst, model_id, ξi, ξ = xi_state

    t0 = time() # to track xi processing duration
        
    # --------------------  TEMP CACHE  --------------------  
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    cache_state = (xi_state, hash(params))
    cached_data = load_cached(cache_state)
    !isnothing(cached_data) && return cached_data

    # Print hello in worker 1
    print_action(xi_state, "STARTING XI")
    
    # --------------------  PARAMS  --------------------  
    ξs::Vector{Float64} = params["ξs"]
    βs::Vector{Float64} = params["βs"]
    ep_alpha::Float64 = params["ep_alpha"]
    ep_epsconv::Float64 = params["ep_epsconv"]
    ep_maxiter::Int = floor(Int, params["ep_maxiter"])
    ep_epoch = params["ep_epoch"]
    stoi_err_th = params["stoi_err_th"]


    # --------------------  PREPARING MODEL  --------------------  
    model::MetNet = deepcopy(ec_models[model_id])
    m, n = size(model)
    obj_idx = rxnindex(model, OBJ_IDER)
    
    # intake info
    intake_info = HG.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossios paper, (see README)
    apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)

    # --------------------  FBA  --------------------  
    fbaout = fba(model, obj_idx)

    # storing
    xi_data = Dict() # local store
    xi_data[(ξ, :net)] = compress_model(model)
    xi_data[(ξ, :fba)] = fbaout
    
    
    # --------------------  BETA LOOP  --------------------  
    # THIS CAN TAKE A WHILE!!!
    βv = zeros(n)
    show_t = time()
    beta_seed = nothing
    for (βi, β) in βs |> enumerate
        
        # --------------------  TEMP CACHE  --------------------  
        # I do not check any cache consistency, so delete the temporal caches if you
        # don't trust them
        # current state
        beta_cache_state = (xi_state, β, hash(params))
        cached = load_cached(beta_cache_state)
        epout = isnothing(cached) ? beta_seed : cached

        best_beta_cache_state = (beta_cache_state, "BEST")
        dat = load_cached(best_beta_cache_state)
        best_epout, best_err = isnothing(dat) ? (nothing, Inf) : dat
        
        # --------------------  TRY MAXENT-EP  --------------------  
        try
            
            # --------------------  MAXENT-EP WHILE LOOP  --------------------  
            # I break down ep in several epoch to be able to cache partial results
            βv[obj_idx] = β
            curr_iter = 0
            while (curr_iter == 0 || isnothing(epout)) || # First time
                (ep_maxiter > curr_iter && # Till maxiter 
                    epout.status != :converged) # Or converged
                
                epout = testing ? empty_epout(model) :
                    maxent_ep(model, 
                            alpha = ep_alpha, 
                            beta_vec = βv, 
                            epsconv = ep_epsconv, 
                            maxiter = ep_epoch,
                            solution = epout, 
                            verbose = false)
                
                curr_iter += epout.iter

                # stoi error
                abs_norm_err = norm_abs_stoi_err(model, epout) .|> abs
                max_err = abs_norm_err |> maximum
                min_err = abs_norm_err |> minimum
                mean_err = abs_norm_err |> mean
                
                # backup best solution 
                if max_err < min(best_err, stoi_err_th)
                    best_epout = epout
                    best_err = max_err
                    save_cache((best_epout, max_err), best_beta_cache_state)
                end

                # Show progress # TODO epout.iter > 10, make this settable
                if βi == 1 || βi == length(βs) || βi % upfrec == 0 || 
                        epout.iter > 10 || epout.status == :converged ||
                        (time() - show_t) > 30 # update if more than 30s
                    show_t = time()
                    
                    exp_av = Rd.val(:μ, stst)
                    fba_av = av(model, fbaout, obj_idx)
                    ep_av = av(model, epout, obj_idx)
                    ep_status = epout.status
                    ep_iter = curr_iter

                    elapsed = time() - t0

                    print_progress(stst, ξi, ξs, ξ, βi, βs, β, exp_av, fba_av, ep_av, 
                        ep_alpha, ep_epsconv, ep_maxiter, ep_status, ep_iter, 
                        min_err, mean_err, max_err,
                        elapsed)
                end
                
                # caching
                save_cache(epout, beta_cache_state)
                
            end # EP While loop

        catch err

            err isa InterruptException && rethrow(err)
            
            # using backup
            isnothing(best_epout) && rethrow(err)
            print_action(beta_cache_state, "ERROR DURING EP", string_err(err))
            epout = best_epout
            
        end # try

        # save for the next beta
        beta_seed = epout

        # storing
        xi_data[(ξ, β, :ep)] = epout

        # caching
        save_cache(epout, beta_cache_state)

    end # β loop
    
    # --------------------  FINISHING --------------------   
    save_cache(xi_data, xi_state)
    
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

## Processing
println("\n\n ------------------ PROCESSING ------------------\n\n")
for (model_id, ec_model) in ec_models
    
    println("\n\n ------------------ $model_id ------------------\n")

    # loading cache
    state = (model_id, hash(params))
    cached = load_cached(state)
    !isnothing(cached) && continue
    
    # this can take a while!!!
    # A, B, C steady states have the same initial conditions
    # ststs_ = [stst for stst in Rd.ststs if stst != "B" && stst != "C"]
    ststs_ = testing ? Rd.ststs[1:1] : Rd.ststs
    
    data_dict = Dict()
    println("Ststs: ", ststs_)
    pmap(ststs_) do stst
        exp_state = (model_id, stst)
        data_dict[stst] = process_exp(exp_state)
    end
    
    save_cache(data_dict, state)
end

## Saving
println("\n\n ------------------ SAVING ------------------\n\n")
boundles = Dict()
for (model_id, ec_model) in ec_models

    # loading cache
    state = (model_id, hash(params))
    cached = load_cached(state)
    isnothing(cached) && continue

    boundles[model_id] = cached
end
file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
tagsave(file, Dict(DATA_KEY => boundles))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")

## Delete Temp Caches
# Do not forget to run this if you change any parameter
println("\n\n ------------------ DELETING CACHE ------------------\n\n")
# delete_temp_caches() # Test
println()
flush(stdout)