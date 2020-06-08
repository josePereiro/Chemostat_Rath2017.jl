# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
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
# using Plots
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

@everywhere notebook_name = "fva_pp_base_model_maxent_ep_v1";

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
            exp_av, fba_av, ep_av, ep_alpha, ep_epsconv, ep_maxiter,
            elapsed)
    Core.println("Worker: $wid stst: $stst progress at $(Time(now())) ----------------------------")
    Core.println("\txi: [$ξi/ $(length(ξs))] beta: [$βi/ $(length(βs))]")
    Core.println("\t  ----------------- --------------")
    Core.println("\tmodel xi:           $ξ")
    Core.println("\tmodel beta:         $β")
    Core.println("\tmodel ep_alpha:     $ep_alpha")
    Core.println("\tmodel ep_epsconv:   $ep_epsconv")
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

@everywhere temp_cache_file(stst) = 
    joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___temp_cache___stst_$(stst).jls")

@everywhere temp_cache_file(stst, xi::Int) = 
    joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___temp_cache___stst_$(stst)_xi_$(xi).jls")

# ---
# ## Params
# ---

# +
@everywhere begin
    
    params = Dict()
    
    # model iders
    params["obj_ider"] = "BIOMASS"
    params["cost_ider"] = "enzyme_solvent_capacity"

    # ep
    params["ep_alpha"] = 1e9
    params["ep_epsconv"] = 1e-5
    params["ep_maxiter"] = 5e3

end
# -

# ---
# ## work functions
# ---

# +
@everywhere function process_stst(stst, upfrec = 10)
    
    # I will cache temporally the results of this function 
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    tcache_file = temp_cache_file(stst)

    # If cached return 
    if isfile(tcache_file)
        
        (stst, boundle) =  deserialize(tcache_file)
        
        # Print in worker 1
        remotecall_wait(print_return_cached, 1, myid(), stst, tcache_file)
        
        return (stst, boundle)
    end
    
    # prepare params
    # Change here how many xi to model, you should always include the experimental xis
    ξs = range(10, 210, length = 3)
    ξs = [ξs; Rd.val("ξ", Rd.ststs)] |> collect |> unique |> sort |> reverse
    
    
    # Change here how many betas to model
    # The beta range is set up by trial and error
    βs = range(10.0^5, 10.0^5.4, length = 30) |> collect |> unique |> sort |> reverse

    
    # Print hello in worker 1
    remotecall_wait(print_stst_hello, 1, myid(), stst, ξs, βs)
    
    boundle = Ch.Utils.ChstatBoundle()
    
    for (ξi, ξ) in ξs |> enumerate 
        
        xi_data = process_xi(stst, ξi, ξs, βs, upfrec)
        
        boundle_xi_data!(boundle, ξ, βs, xi_data)
        
    end
    
    # Catching 
    serialize(tcache_file, ((stst, boundle)))
    
    # Printing in 1
    remotecall_wait(print_stst_good_bye, 1, myid(), stst, tcache_file)
    
    return (stst, boundle)
    
end

# +
@everywhere function process_xi(stst, ξi, ξs, βs, upfrec)
    
    ξ = ξs[ξi]
        
    # If cached load 
    xi_tcache_file = temp_cache_file(stst, ξi)
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
    ep_epsconv = params["ep_epsconv"]
    ep_maxiter = floor(Int, params["ep_maxiter"])

    obj_ider = params["obj_ider"]
    cost_ider = params["cost_ider"]
    intake_info = M.stst_base_intake_info(stst)


    # prepare model
    model = deserialize(M.FVA_PP_BASE_MODEL_FILE);
    obj_idx = Ch.Utils.rxnindex(model, obj_ider)
    

    # Chemostat steady state constraint, see Cossios paper, (see README)
    Ch.SteadyState.apply_bound!(model, ξ, intake_info)

    # atat demand
    # The working cell line is a production line
    M.add_a1at_synthesis!(model, stst)

    # maxent-fba
    fbaout = Ch.FBA.fba(model, obj_ider, cost_ider)

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

            seed_epout = epout isa Ch.Utils.AbstractOut ? epout : nothing

            # storing
            xi_data[(ξ, β, :ep)] = epout

            # Info
            ep_av = Ch.Utils.av(model, epout, obj_idx)

            # Print progress in worker 1
            show_progress = βi == 1 || βi == length(βs) || βi % upfrec == 0
            show_progress && 
                remotecall_wait(print_progress, 1, myid(), stst, 
                    ξi, ξs, ξ, βi, βs, β, 
                    exp_av, fba_av, ep_av, 
                    ep_alpha, ep_epsconv, ep_maxiter,
                    time() - t0)
            
        catch err
            # Print error in 1
            remotecall_wait(print_error, 1, myid(), stst, ξi, βi, err)
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

    Ch.Utils.add_data!(boundle, ξ, :net, xi_data[(ξ, :net)])
    Ch.Utils.add_data!(boundle, ξ, :fba, xi_data[(ξ, :fba)])
    
    for β in βs
        Ch.Utils.add_data!(boundle, ξ, β, :ep, xi_data[(ξ, β, :ep)])
    end
end

# ### Parallel loop

# this can take a while!!!
# A, B, C steady states have the same initial conditions
ststs_ = [stst for stst in Rd.ststs if stst != "B" && stst != "C"]
println("Ststs: ", ststs_)
remote_results = pmap(process_stst, ststs_);

# ### Delete Temp Caches
# It is safer to run this last, but I have not enought memory to have redundant results

# +
# Do not forget to run this if you change any parameter
for stst in Rd.ststs
    
    # Steady state cache
    tcache_file = temp_cache_file(stst)
    if isfile(tcache_file)
        rm(tcache_file, force = true)
        println(relpath(tcache_file), " deleted!!!"); flush(stdout)
    end
    
    # Xi caches
    for xi in 1:10_000 # More than this, I dont think so
        tcache_file = temp_cache_file(stst, xi)
        if isfile(tcache_file)
            rm(tcache_file, force = true)
            println(relpath(tcache_file), " deleted!!!"); flush(stdout)
        end
    end
    
end
# -

## Saving
file_ = joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___boundles.jls")
serialize(file_, (params, remote_results))
println(relpath(file_), " created!!!")

# +
# stst, boundle = remote_results[1];

# +
# # using Plots
# p = Plots.plot()
# Ch.Plots.plot_marginal!(p, boundle, 1, 1, [:ep, :fba], "BIOMASS")
# -









# +
# for stst in Rd.ststs
    
#     ξs = range(10, 250, length = 3)
#     ξs = [ξs; Rd.val("ξ", Rd.ststs)] |> unique |> sort

    
#     # The beta range is set up by trial and error
#     βs = floor.(Ch.Utils.logspace(5, 5.6, 30)) |> sort
    
#     # both have the same initial conditions than C
# #     stst == "A" && continue
# #     stst == "B" && continue
    
#     # rath steady state
#     meta_info["stst"] = stst
    
#     # the maximum beta at which we get results
#     max_feasible_β = nothing 
    
#     # Taken possible cached boundle
#     boundle = haskey(boundles, stst) ? boundles[stst] : Ch.Utils.ChstatBoundle()

#     intake_info = M.stst_base_intake_info(stst)

#     # Computing
#     ξs = [Rd.val("ξ", stst)]
#     ξs_str = round.(ξs, digits = 3)
#     for (ξi, (ξ, ξstr)) in enumerate(zip(ξs, ξs_str))        

#         # prepare model
#         model = deserialize(M.FVA_PP_BASE_MODEL_FILE);
#         obj_idx = Ch.Utils.rxnindex(model, obj_ider)

#         # Chemostat steady state constraint, see Cossios paper, (see README)
#         Ch.SteadyState.apply_bound!(model, ξ, intake_info)
        
#         # atat demand
#         M.add_a1at_synthesis!(model, stst)

#         # maxent-fba
#         fbaout = Ch.FBA.fba(model, obj_idx)

#         # boundling
#         Ch.Utils.add_data!(boundle, ξ, :fba, fbaout)
#         Ch.Utils.add_data!(boundle, ξ, :net, model)
        
#         # seed solution
#         seed_epout = nothing

#         βv = zeros(size(model, 2))
#         β_iter_ = zip(sort(βs, rev = true), sort(βs_str, rev = true))
#         for (βi, (β, βstr)) in β_iter_ |> enumerate
            
#             # avoid error runs
#             !isnothing(max_feasible_β) && β > max_feasible_β && continue
            
#             # avoid redo
#             haskey(boundle, ξ, β, :ep) && continue
            
#             # Info
#             seed_str = isnothing(seed_epout) ? "No seed" : "Seed mu: $(Ch.Utils.av(model, seed_epout, obj_ider))"
#             println("Doing stst: $stst, xi [$ξi / $(length(ξs))]: $ξstr  beta [$βi / $(length(βs))]: $βstr    $(seed_str)"); flush(stdout)

#             try
#                 # maxent-ep
#                 βv[obj_idx] = β
#                 @time epout = Ch.MaxEntEP.maxent_ep(model, α = α, βv = βv, 
#                     epsconv = epsconv_, solution = seed_epout, verbose = true); flush(stdout)
                
#                 # boundling
#                 Ch.Utils.add_data!(boundle, ξ, β, :ep, epout)

#                 # Info
#                 fba_mu = Ch.Utils.av(boundle, ξ, :fba, obj_ider)
#                 ep_mu = Ch.Utils.av(boundle, ξ, β, :ep, obj_ider)
#                 println("\tfba mu: $(fba_mu)\tep_mu: $(ep_mu)"); flush(stdout)
                
#                 # seed and max_feasible_β
#                 if isnothing(max_feasible_β) && epout isa Ch.Utils.EPout
#                     max_feasible_β = β
#                     println("Max feasible beta[$βi]: $β")
#                 end
#                 seed_epout = epout
                
#             catch err
#                 err isa InterruptException && rethrow(err)
                
#                 # boundling
#                 Ch.Utils.add_data!(boundle, ξ, β, :ep, err)
#                 println("ERROR Doing xi: $ξstr  beta: $βstr, $(err)"); flush(stdout)
#                 continue
#             end # try

#         end # β loop
#         println()
        
#     end # ξ loop
#     println()
    
#     boundles[stst] = boundle
#     # Catching
#     cache_file = joinpath(M.MODEL_CACHE_DATA_DIR, "$(notebook_name)___boundle___$(stst).jls");
#     serialize(cache_file, (meta_info, boundle))
#     @assert isfile(cache_file)
#     println("created $(relpath(cache_file)) !!!")

# end # stst loop
# println("Finitooooooooooooo!!! ", " "^50)

# +
# ξs = copy(boundle.ξs)
# ξs_str = round.(ξs, digits = 3);
# println("ξs ($(length(ξs))) from ", minimum(ξs), " to ", maximum(ξs))

# βs = copy(boundle.βs);
# βs_str = round.(βs, digits = 3);
# println("βs ($(length(βs))) from ", minimum(βs), " to ", maximum(βs))
# -

# #### Biomass checking

# +
# ider = obj_ider
# p = Plots.plot(title = ider, 
#     xlabel = "xi", xscale = :linear,
#     ylabel = "mu", yaxis = [-0.01, 0.08])
# avs = Ch.Utils.av(boundle, ξs, :fba, ider)
# Plots.scatter!(boundle.ξs, avs,  label = "fba", lw = 3, ms = 5)
# Plots.scatter!(Rd.val("ξ", Rd.ststs), Rd.val("μ", Rd.ststs), label = "exp", ms = 5, 
#     xerr = Rd.err("ξ", Rd.ststs), yerr = Rd.err("μ", Rd.ststs), color = :black)
# -

# ---
# ## CheckPlots
# --- 

# +
# ider = obj_ider
# eps_ = minimum(βs) == 0 ? 0.0 : 0.1
# p = Plots.plot(title = ider, legend = :topleft, 
#     xlabel = "beta",
#     ylabel = "flx",
#     xscale = :log10)
# Plots.scatter!(p, [], [], label = "ep", color = :black) # legend
# for (βstr, β) in zip(βs_str, βs)
#     ep_av = Ch.Utils.av(boundle, ξs[1], β, :ep, ider)
#     Plots.scatter!(p, [β + eps_], [ep_av], ms = 2, label = "", color = :black)   
# end
# fba_av = Ch.Utils.av(boundle, ξs[1], :fba, ider)
# Plots.plot!(p, x -> fba_av, minimum(βs) + eps_, maximum(βs),
#     lw = 3, ls = :dot, color = :red,
#     label = "fba")
# Plots.plot!(p, x -> Rd.val("μ", stst), minimum(βs) + eps_, maximum(βs),
#     lw = 3, ls = :dash, label = "exp E")
# p
