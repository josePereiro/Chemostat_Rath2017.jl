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

# ---
# ## Description

# This script produce the plots from the date generated from [fva_pp_tINIT_models_maxent_ep](./fva_pp_tINIT_models_maxent_ep.jl).

# +
using DataFrames
using Serialization
using Dates
using StatsBase

using Plots
pyplot();
# -

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
HG = Chemostat_Rath2017.HumanGEM
tIG = Chemostat_Rath2017.tINIT_GEMs;

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# ## MaxEnt EP results

# +
#TODO: outhomatize this script for all the networks
# -

dat_files = filter(readdir(tIG.MODEL_PROCESSED_DATA_DIR)) do file
    startswith(file, "fva_pp_tINIT_models_maxent_ep___")
end
dat_files = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, dat_files);

dat = deserialize(dat_files[1])
dat_id = dat.id
boundles = Dict()
for dat in dat.res
    stst, boundle = dat
    boundles[stst] = boundle
end

obj_ider = "biomass_human"
good_βs = Dict()
for (stst, boundle) in boundles
    good_βs[stst] = filter(boundle.βs) do β
        dat = boundle[1, β, :ep]
        dat isa Ch.Utils.EPout && dat.status == :converged
    end
end
println("good betas: ", length(good_βs))

println("\nExp_beta")
exp_βs = Dict()
for (stst, boundle) in boundles
    exp_μ = Rd.val(:μ, stst)
    βs = good_βs[stst]
    exp_β = βs |> first
    for β in βs
        ep_μ = Ch.Utils.av(boundle, 1, β, :ep, obj_ider)
        last_ep_μ = Ch.Utils.av(boundle, 1, exp_β, :ep, obj_ider)
        if abs(ep_μ - exp_μ) < abs(last_ep_μ - exp_μ)
            exp_β = β
        end
    end
    ep_μ = Ch.Utils.av(boundle, 1, exp_β, :ep, obj_ider)
    println("stst: ", stst)
    println("good_betas: ", length(βs))
    println("exp beta: ", exp_β)
    println("exp μ: ", exp_μ)
    println("ep μ:  ", ep_μ)
    println()
    exp_βs[stst] = exp_β
end

# ---
# ## Plots
# ---

colors = [:red, :yellow, :blue, :green, :orange, :black];

## Stoi err
p = plot(title = dat_id, xlabel = "beta", ylabel = "abs norm stoi error")
for (i, (stst, boundle)) in boundles |> enumerate
    βs = good_βs[stst] |> sort
    metnet = boundle[1, :net]
    epouts = boundle[1, βs, :ep]
    errs = map(epouts) do epout 
        Ch.Utils.norm_abs_stoi_err(metnet, epout)
    end
    exp_β = exp_βs[stst]
    vline!([exp_β], ls = :dot, lw = 2, color = colors[i], label = "")
    plot!(βs, mean.(errs), label = "", color = colors[i], lw = 2, ls = :dash)
    plot!(βs, maximum.(errs), label = "", color = colors[i], lw = 2)
    scatter!([],[], label = stst, marker = :square, ms = 30, color = colors[i])
end
plot!([],[], ls = :dash, color = :black, label = "mean")
plot!([],[], ls = :dot, color = :black, label = "exp_beta")
plot!([],[], color = :black, label = "max")

file = joinpath(tIG.MODEL_FIGURES_DATA_DIR, "fva_pp___$(dat_id)___norm_stoi_err.png")
# savefig(p, file);
println(relpath(file), " created!!!")

# correlations at exp_beta
fluxs = [obj_ider; [HG.exch_met_map[HG.mets_map[rath_met]] for rath_met in Rd.msd_mets]];
ep_p = plot(title = "$dat_id\nEP", xlabel = "exp", ylabel = "ep", 
    ylim = [-5, 25]
)
ep_μ_p = plot(title = "$dat_id\nμ EP", xlabel = "exp", ylabel = "ep")
fba_p = plot(title = "$dat_id\nFBA", xlabel = "exp", ylabel = "fba", 
    ylim = [-0.2, 0.2]
)
for (i, (stst, boundle)) in boundles |> enumerate
    β = exp_βs[stst]
    for rath_id in ["μ"; Rd.msd_mets]
        exp_q = rath_id == "μ" ? 
            Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
        exp_err = rath_id == "μ" ? 
            Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)
        model_id = rath_id == "μ" ? 
            obj_ider : HG.exch_met_map[HG.mets_map[rath_id]]
        ep_q = Ch.Utils.av(boundle, 1, β, :ep, model_id)
        ep_err = sqrt(Ch.Utils.va(boundle, 1, β, :ep, model_id))
        scatter!(ep_p, [exp_q], [ep_q], 
            xerr = [exp_err], yerr = [ep_err],
            label = "", color = colors[i], ms = 10)
        rath_id == "μ" && scatter!(ep_μ_p, [exp_q], [ep_q], 
            xerr = [exp_err], yerr = [ep_err],
            label = "", color = colors[i], ms = 10)
        
        fba_q = Ch.Utils.av(boundle, 1, :fba, model_id)
        fba_err = sqrt(Ch.Utils.va(boundle, 1, :fba, model_id))
        scatter!(fba_p, [exp_q], [fba_q], 
            xerr = [exp_err], 
            label = "", color = colors[i], ms = 10)
                
        # info
#         tol = 10
#         abs(ep_q/exp_q) > tol && println(stst, " ", rath_id, " > $tol difference")
    end
end
plot(ep_μ_p, ep_p, fba_p, size = [1200, 400], layout = (1, 3))

# +
# correlation as function of beta
ep_ps = []


for β in intersect(values(good_βs)...)
    
    βstr = round(β; digits = 3)
    
    ep_p = plot(title = "$dat_id\nEP-beta: $βstr", xlabel = "exp", ylabel = "ep", 
            ylim = [-5, 40]
        )
    for (i, (stst, boundle)) in boundles |> enumerate

        for rath_id in ["μ"; Rd.msd_mets]
            exp_q = rath_id == "μ" ? 
                Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
            exp_err = rath_id == "μ" ? 
                Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)
            model_id = rath_id == "μ" ? 
                obj_ider : HG.exch_met_map[HG.mets_map[rath_id]]
            ep_q = Ch.Utils.av(boundle, 1, β, :ep, model_id)
            ep_err = sqrt(Ch.Utils.va(boundle, 1, β, :ep, model_id))
            scatter!(ep_p, [exp_q], [ep_q], 
                xerr = [exp_err], yerr = [ep_err],
                label = "", color = colors[i], ms = 10)
        end
        
    end
    push!(ep_ps, ep_p)
end
# -

ps_ = ep_ps[1:3:end]
p = plot(ps_..., size = (3000,300), layout = (1, length(ps_)))

#Marginals
function plot_marginals(rath_id; kwargs...)
    ps = []
    for (i, (stst, boundle)) in boundles |> enumerate
        β = exp_βs[stst]

        # averages
        exp_q = rath_id == "μ" ? 
                Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
        exp_err = rath_id == "μ" ? 
            Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)
        model_id = rath_id == "μ" ? 
            obj_ider : HG.exch_met_map[HG.mets_map[rath_id]]
        ep_q = Ch.Utils.av(boundle, 1, β, :ep, model_id)
        ep_err = sqrt(Ch.Utils.va(boundle, 1, β, :ep, model_id))
        
        p = plot(title = model_id, xlabel = "pdf", ylabel = "flux"; kwargs...)
        vline!([ep_q], ls = :dash, label = "", color = :black, lw = 3)
        vline!([exp_q], ls = :dash, label = "", color = :green, lw = 3)
        Ch.Plots.plot_marginal!(p, boundle, 1, β, [:fba, :ep], model_id; label = "")

        push!(ps, p)
    end
    return plot(ps..., size = [1500,300], layout = grid(1, length(ps)))
end

# plot_marginals("μ", xlim = [0.0, 0.05])
# plot_marginals("GLN", xlim = [-0.5, 1])
plot_marginals("GLC", xlim = [-0.5, 4])

import LinearAlgebra

base_model = deserialize(HG.BASE_MODEL_FILE);

#=
    So what do you do when you have an ill-conditioned matrix? 
    You reformulate your problem by doing something called preconditioning. 
=#
LinearAlgebra.cond(base_model.S) # = 2.8150143905656065e20


