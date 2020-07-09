# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl:light,ipynb
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

# ### Precompaling in master worker first

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
M = Chemostat_Rath2017.MODEL1105100000;

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# ### Loading everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES && addprocs(NO_CORES - 1)
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
M = Chemostat_Rath2017.MODEL1105100000;
    
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env()
    
end
# -

# ---
# ## META
# ---

@everywhere notebook_name = "fva_pp_base_model_maxent_ep_epsconv_study";

# ---
# ## Print functions
# ---

@everywhere function print_stst_hello(wid, stst, ξs, βs)
    Core.println("Worker $wid starting stst $stst at $(Time(now())) ----------------------------")
    Core.println("\txis:   $ξs")
    Core.println("\tbetas: $βs")
    Core.println()
    flush(stdout);
end

@everywhere function print_xi_hello(wid, stst, ξ, βs)
    Core.println("Worker $wid stst $stst starting xi $ξ at $(Time(now())) ----------------------------")
    Core.println("\tbetas: $βs")
    Core.println()
    flush(stdout);
end

@everywhere function print_progress(wid, stst, 
            ξi, ξs, ξ,  βi, βs, β, 
            exp_av, fba_av, 
            ep_av, ep_alpha, ep_epsconv, ep_maxiter, ep_status, ep_iter,
            elapsed)
    Core.println("Worker: $wid stst: $stst progress at $(Time(now())) ----------------------------")
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

@everywhere function print_stst_good_bye(wid, stst, tcache_file)
    Core.println("Worker $wid finished stst $stst at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

@everywhere function print_xi_good_bye(wid, stst, xi, tcache_file)
    Core.println("Worker $wid finished stst $stst xi $xi at $(Time(now())) ----------------------------")
    Core.println("\tXi temp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

# +
@everywhere function print_return_cached(wid, stst, tcache_file)
    Core.println("Worker $wid returns cached stst $stst at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

@everywhere function print_return_cached(wid, stst, xi, tcache_file)
    Core.println("Worker $wid stst $stst load cached xi $xi at $(Time(now())) ----------------------------")
    Core.println("\tXi temp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end
# -

@everywhere function print_error(wid, stst, xi, beta, err)
    Core.println("Worker $wid stst $stst xi $xi beta $beta ERROR at $(Time(now())) ----------------------------")
    Core.println("\t", err)
    Core.println()
    flush(stdout);
end

# ### temp cache file

@everywhere temp_cache_file(state...) = 
    joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___temp_cache___state_$(hash(state)).jls")

# ---
# ## Params
# ---

# +
@everywhere begin
    
    params = Dict()
    params["stst"] = "A"
    
    # model iders
    params["obj_ider"] = "BIOMASS"
    params["cost_ider"] = "enzyme_solvent_capacity"
    params["cost_met_ider"] = "enzyme_c"
    

    # ep
    params["ep_alpha"] = 1e9
    params["ep_epsconvs"] = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    params["ep_maxiter"] = 1e4
    # beta values that force the approximation of the current model
    # to the experimental objective (μ). 
    # The values come from previows experiments
    params["closest_βs"] = Dict(
        "B"=>1.57347e5,
        "A"=>1.67774e5,
        "C"=>1.41707e5,
        "F01"=>1.67774e5,
        "D"=>1.67774e5,
        "E"=>1.62561e5)
        
end
# -

# ---
# ## work functions
# ---

# +
@everywhere function process_state(state, upfrec = 10)
    
    # extract state
    stst, epsconv = state
    
    # I will cache temporally the results of this function 
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    tcache_file = temp_cache_file(state)

    # If cached return 
    if isfile(tcache_file)
        
        (state, boundle) =  deserialize(tcache_file)
        
        # Print in worker 1
        remotecall_wait(print_return_cached, 1, myid(), stst, tcache_file)
        
        return (state, boundle)
    end
    
    # prepare params
    # Change here how many xi to model, you should always include the experimental xis
    ξs = [Rd.val("ξ", stst)] 
    
    # Change here how many betas to model
    # The beta range is set up by trial and error
    closest_β = params["closest_βs"][stst]
    βs = [closest_β]
    
    # Print hello in worker 1
    remotecall_wait(print_stst_hello, 1, myid(), stst, ξs, βs)
    
    # xi_data
    xis_data = map((ξi) -> process_xi(state, ξi, ξs, βs, upfrec), eachindex(ξs))
    
    # Boundling
    boundle = Ch.Utils.ChstatBoundle()
    foreach(eachindex(ξs)) do ξi
        boundle_xi_data!(boundle, ξs[ξi], βs, xis_data[ξi])
    end

    # Catching 
    serialize(tcache_file, (state, boundle))
    
    # Printing in 1
    remotecall_wait(print_stst_good_bye, 1, myid(), stst, tcache_file)
    
    return (state, boundle)
    
end

# +
@everywhere function process_xi(state, ξi, ξs, βs, upfrec)
    
    # extract state
    stst, ep_epsconv = state
    
    ξ = ξs[ξi]
        
    # If cached load 
    xi_tcache_file = temp_cache_file(state, ξ)
    if isfile(xi_tcache_file)
        xi_data =  deserialize(xi_tcache_file)

        # Print in worker 1
        remotecall_wait(print_return_cached, 1, myid(), stst, ξi, xi_tcache_file)

        return xi_data
    end
    
    # Info
    # Print hello in worker 1
    remotecall_wait(print_xi_hello, 1, myid(), stst, ξi, βs)
    
    t0 = time() # to track xi processing duration

    # temporally store
    xi_data = Dict()

    # Params
    ep_alpha = params["ep_alpha"]
    ep_maxiter = floor(Int, params["ep_maxiter"])

    obj_ider = params["obj_ider"]
    cost_ider = params["cost_ider"]
    cost_met_ider = params["cost_met_ider"]
    intake_info = M.stst_base_intake_info(stst)


    # prepare model
    model = deserialize(M.FVA_PP_BASE_MODEL_FILE);
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)
    
    # delete cost
    cost_met_idx = Ch.Utils.metindex(model, cost_met_ider)
    model.S[cost_met_idx, :] .= 0.0

    # Chemostat steady state constraint, see Cossios paper, (see README)
    Ch.SteadyState.apply_bound!(model, ξ, intake_info)

    # atat demand
    # The working cell line is a production line
    M.add_a1at_synthesis!(model, stst)

    # maxent-fba
    fbaout = Ch.LP.fba(model, obj_ider, cost_ider)
    

    # storing
    xi_data[(ξ, :net)] = model
    xi_data[(ξ, :fba)] = fbaout
    
    # Info
    exp_av = Rd.val(:μ, stst)
    fba_av = Ch.Utils.av(model, fbaout, obj_idx)
    
    βv = zeros(size(model, 2))
    seed_epout = nothing
    for (βi, β) in βs |> enumerate
        
        try
            # epout
            βv[obj_idx] = β
            epout = Ch.MaxEntEP.maxent_ep(model, 
                alpha = ep_alpha, 
                beta_vec = βv, 
                epsconv = ep_epsconv, 
                maxiter = ep_maxiter,
                solution = seed_epout, 
                verbose = false)

            seed_epout = epout

            # storing
            xi_data[(ξ, β, :ep)] = epout

            # Info
            ep_av = Ch.Utils.av(model, epout, obj_idx)

            # Print progress in worker 1
            show_progress = βi == 1 || βi == length(βs) || βi % upfrec == 0 || epout.iter > 10
            show_progress && 
                remotecall_wait(print_progress, 1, myid(), stst, 
                    ξi, ξs, ξ, βi, βs, β, 
                    exp_av, fba_av, ep_av, 
                    ep_alpha, ep_epsconv, ep_maxiter, epout.status, epout.iter,
                    time() - t0)

        catch err
            # Print error in 1
            remotecall_wait(print_error, 1, myid(), stst, ξi, βi, err)
            # storing
#             xi_data[(ξ, β, :ep)] = err
            rethrow(err)
            err isa InterruptException && rethrow(err)
        end
        
    end
    
    # Catching 
    serialize(xi_tcache_file, xi_data)
    
    # Printing in 1
    remotecall_wait(print_xi_good_bye, 1, myid(), stst, ξi, xi_tcache_file)
    
    return xi_data
    
end
# -

@everywhere function boundle_xi_data!(boundle, ξ, βs, xi_data)

    boundle[ξ, :net] = xi_data[(ξ, :net)]
    boundle[ξ, :fba] = xi_data[(ξ, :fba)]
    
    for β in βs
        boundle[ξ, β, :ep] = xi_data[(ξ, β, :ep)]
    end
end

# ### Processing

# +
# this can take a while!!!
# A, B, C steady states have the same initial conditions
# -

states = Iterators.product(Rd.ststs, params["ep_epsconvs"]) |> collect |> vec;

println()
println("States")
println.(states);
println()
remote_results = pmap(process_state, states);

# ### Delete Temp Caches
# It is safer to run this last, but I have not enought memory to have redundant results

# Do not forget to run this if you change any parameter
cache_dir = dirname(temp_cache_file("bla"))
temp_file_hint = "$(notebook_name)___temp_cache"
for file in readdir(M.MODEL_CACHE_DATA_DIR)
    if startswith(file, temp_file_hint)
        tcache_file = joinpath(cache_dir, file)
        rm(tcache_file, force = true)
        println(relpath(tcache_file), " deleted!!!")
    end
end

# ---
# ## Saving
# ---

file_ = joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___boundles.jls")
serialize(file_, (params, remote_results))
println(relpath(file_), " created!!!")
