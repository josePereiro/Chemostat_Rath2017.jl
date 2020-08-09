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

# This script produce plots from the exctracted data in [maxent_ep___extract_data](./maxent_ep___extract_data.jl).

using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using DataFrames
using Serialization
using Dates
using StatsBase
using JSON
using Plots
pyplot();

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
const ChR = Chemostat_Rath2017
const Rd = ChR.RathData
const HG = ChR.HumanGEM
const tIG = ChR.tINIT_GEMs;

# ---
# ## Load data
# ---

# To see `models_dat` layout see [maxent_ep___extract_data](./maxent_ep___extract_data.jl).

# maxent_ep data files
# TODO: package this
preffix = "fva_pp_tINIT_models_maxent_ep___"
suffix = "___boundles___extracted_data.json"
# TODO: use DrWatson
models_ids = ["GTEx-brain-1", "TCGA-GBM NT-1", "TCGA-GBM TP-1", 
    "TCGA-GBM TR-1", "TCGA-LGG TP-1", "TCGA-LGG TR-1"]
models_dat = Dict()
foreach(readdir(tIG.MODEL_PROCESSED_DATA_DIR)) do file
    if startswith(file, preffix) && endswith(file, suffix)
        id = file[length(preffix) + 1:end - length(suffix)]
        models_dat[id] = JSON.parsefile(joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, file));
    end
end
println("Found: ")
models_dat |> keys  .|> basename .|>  println;

# ---
# ## Plots
# ---

colors = [:red, :yellow, :blue, :green, :orange, :black];

# ## growth vs beta

function plot_growth_vs_beta(dat_id, rath_id = "μ"; kwargs...)
    dat = models_dat[dat_id]
    p = plot(title = "$dat_id\n$rath_id", xlabel = "beta", ylabel = rath_id)
    for (i, stst) in Rd.ststs |> enumerate
        βs = dat[stst]["βs"] |> sort
        ξ = dat[stst]["ξs"] |> first
        
        for β in βs
            ep_dat = dat[stst]["""($ξ, $β, :ep, "$rath_id")"""]
            ep_q = ep_dat["av"] 
            ep_err = sqrt(ep_dat["va"])
            scatter!(p, [β], [ep_q], 
                yerr = [ep_err], 
                label = "", color = colors[i], ms = 10)
        end
        scatter!(p, [],[], label = stst, marker = :square, ms = 30, color = colors[i])
        
    end
    return plot!(p; kwargs...)
end

ps = []
for model_id in models_ids
    ps_ = []
    for rath_id in ["μ"; Rd.msd_mets]
        p = plot_growth_vs_beta(models_ids[2], rath_id)
        push!(ps_, p)
    end
    p = plot(ps_..., size = (length(ps_) * 300, 300), layout = grid(1, length(ps_)))
    push!(ps, p)
end
p = plot(ps...,  size = (2500, length(ps_) * 300), layout = grid(length(ps), 1))

# ## Stoi err

function plot_norm_stoi_err(dat_id; kwargs...)
    dat = models_dat[dat_id]
    p = plot(title = dat_id, xlabel = "beta", ylabel = "abs norm stoi error")
    for (i, stst) in Rd.ststs |> enumerate
        βs = dat[stst]["βs"] |> sort
        exp_β = dat[stst]["exp_β"]
        vline!([exp_β], ls = :dot, lw = 2, color = colors[i], label = "")
        ξ = dat[stst]["ξs"] |> first
#         for ξ in dat[stst]["ξs"]
            mean_err = dat[stst]["""($ξ, "mean_ep_norm_stoi_err")"""]
            plot!(βs, mean_err, label = "", color = colors[i], lw = 2, ls = :dash)
            max_err = dat[stst]["""($ξ, "max_ep_norm_stoi_err")"""]
            plot!(βs, max_err, label = "", color = colors[i], lw = 2)
            scatter!(p, [],[], label = stst, marker = :square, ms = 30, color = colors[i])
#         end
    end

    plot!(p, [],[], ls = :dash, color = :black, label = "mean")
    plot!(p, [],[], ls = :dot, color = :black, label = "exp_beta")
    plot!(p, [],[], color = :black, label = "max")
    
    return plot!(p; kwargs...)
end

ps = []
for models_id in models_ids
    p = plot_norm_stoi_err(models_id, legend = models_id == models_ids[end])
    push!(ps, p)
end
p = plot(ps..., size = (500 * length(ps), 400), layout = grid(1, length(ps)))

# +
# file = joinpath(tIG.MODEL_FIGURES_DATA_DIR, "fva_pp___$(dat_id)___norm_stoi_err.png")
# savefig(p, file);
# println(relpath(file), " created!!!")
# -
# ## correlations at exp_beta

# fluxs = [obj_ider; [HG.exch_met_map[HG.mets_map[rath_met]] for rath_met in Rd.msd_mets]];
function plot_exp_β_corrs(dat_id; kwargs...)
    dat = models_dat[dat_id]
    ep_p = plot(title = "$dat_id\nEP", xlabel = "exp", ylabel = "ep", 
        ylim = [-5, 25]
    )
    ep_μ_p = plot(title = "$dat_id\nμ EP", xlabel = "exp", ylabel = "ep")
    fba_p = plot(title = "$dat_id\nFBA", xlabel = "exp", ylabel = "fba", 
        ylim = [-0.2, 0.2]
    )
    for (i, stst) in Rd.ststs |> enumerate
        ider_map = dat[stst]["ider_map"]
        β  = dat[stst]["exp_β"]
        ξ = dat[stst]["ξs"] |> first
        
        for rath_id in ["μ"; Rd.msd_mets]
            
            # Experimental data
            exp_q = rath_id == "μ" ? 
                Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
            exp_err = rath_id == "μ" ? 
                Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)
            
            # model data
            ep_dat = dat[stst]["""($ξ, $β, :ep, "$rath_id")"""]
            ep_q = ep_dat["av"] 
            ep_err = sqrt(ep_dat["va"])
            
            scatter!(ep_p, [exp_q], [ep_q], 
                xerr = [exp_err], yerr = [ep_err],
                label = "", color = colors[i], ms = 10)
            rath_id == "μ" && scatter!(ep_μ_p, [exp_q], [ep_q], 
                xerr = [exp_err], yerr = [ep_err],
                label = "", color = colors[i], ms = 10)

            fba_dat = dat[stst]["""($ξ, :fba, "$rath_id")"""]
            fba_q = fba_dat["av"] 
            fba_err = 0.0
            scatter!(fba_p, [exp_q], [fba_q], 
                xerr = [exp_err], 
                label = "", color = colors[i], ms = 10)

    # #         info
    #         tol = 100
    #         abs(ep_q/exp_q) > tol && println(stst, ": ", rath_id, " > $tol difference")
        end
    end
    p = plot(ep_μ_p, ep_p, fba_p, size = [1200, 400], layout = (1, 3))
    plot!(p; kwargs...)
end

ps = []
for models_id in models_ids
    p = plot_exp_β_corrs(models_id, titlefont = 10)
    push!(ps, p)
end
p = plot(ps..., size = (1000, 300 * length(ps)), layout = grid(length(ps), 1))

# ## correlation as function of beta

function get_ep_βs_corrs(dat_id)
    dat = models_dat[dat_id]
    ep_ps = []
    βs = [dat[stst]["βs"] for stst in Rd.ststs]
    for β in intersect(βs...)
        βstr = round(β; digits = 3)

        ep_p = plot(title = "$dat_id\nEP-beta: $βstr", xlabel = "exp", ylabel = "ep", 
                ylim = [-5, 40]
            )
        for (i, stst) in Rd.ststs |> enumerate
            
            ξ = dat[stst]["ξs"] |> first

            for rath_id in ["μ"; Rd.msd_mets]
                
                # Experimental data
                exp_q = rath_id == "μ" ? 
                    Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
                exp_err = rath_id == "μ" ? 
                    Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)
                                
                # model data
                ep_dat = dat[stst]["""($ξ, $β, :ep, "$rath_id")"""]
                ep_q = ep_dat["av"] 
                ep_err = sqrt(ep_dat["va"])

                scatter!(ep_p, [exp_q], [ep_q], 
                    xerr = [exp_err], yerr = [ep_err],
                    label = "", color = colors[i], ms = 10)
            end

        end
        push!(ep_ps, ep_p)
    end
    
    return ep_ps
end

# ps_ = get_ep_βs_corrs(models_ids[1])
# ps_ = ps_[1:3:end]
# p = plot(ps_..., size = (3000,300), layout = (1, length(ps_)))
ps = []
for models_id in models_ids
    ps_ = get_ep_βs_corrs(models_id)
    ps_ = ps_[1:3:end]
    p = plot(ps_..., size = (3000,300), layout = (1, length(ps_)))
    push!(ps, p)
end
p = plot(ps..., size = (3000, 300 * length(ps)), layout = grid(length(ps), 1))

#Marginals
function plot_marginals(dat_id, rath_id; kwargs...)
    dat = models_dat[dat_id]
    
    ps = []
    for (i, stst) in Rd.ststs |> enumerate
        β = dat[stst]["exp_β"]
        ξ = dat[stst]["ξs"] |> first
        
        # Experimental data
        exp_q = rath_id == "μ" ? 
            Rd.val("μ", stst) : -Rd.val("q" * rath_id, stst)
        exp_err = rath_id == "μ" ? 
            Rd.err("μ", stst) : Rd.err("q" * rath_id, stst)

        # model data
        ep_dat = dat[stst]["""($ξ, $β, :ep, "$rath_id")"""]
        ep_q = ep_dat["av"] 
        ep_err = sqrt(ep_dat["va"])
        lb, ub = dat[stst]["""($ξ, :bounds, "$rath_id")"""]
        ep_μ = ep_dat["μ"]
        ep_σ = ep_dat["σ"]

        p = plot(;kwargs...)
        Ch.Plots.plot_marginal!(p, ep_μ, ep_σ, lb, ub)
        vline!([ep_q], ls = :dash, label = "", color = :black, lw = 3)
        vline!([exp_q], ls = :dash, label = "", color = :green, lw = 3)
        p = plot!(p, title = "$dat_id\n$rath_id ($stst)", xlabel = "", ylabel = "")
        push!(ps, p)
    end
    return plot(ps..., size = [1500,300], layout = grid(1, length(ps)))
end

plot_marginals(models_ids[1], "μ", xlim = [0.0, 0.05])

ps = []
push!(ps, plot_marginals("μ", xlim = [0.0, 0.05]))
push!(ps, plot_marginals("GLC", xlim = [-0.5, 10]))
push!(ps, plot_marginals("LAC", xlim = [-0.5, 35]))
push!(ps, plot_marginals("GLN", xlim = [-0.5, 1.5]))
push!(ps, plot_marginals("NH4", xlim = [-0.5, 10]))
push!(ps, plot_marginals("GAL", xlim = [-0.5, 15]))
push!(ps, plot_marginals("PYR", xlim = [-0.5, 1000]))
push!(ps, plot_marginals("GLU", xlim = [-0.5, 13]))
push!(ps, plot_marginals("ALA", xlim = [-0.5, 13]))
push!(ps, plot_marginals("ASP", xlim = [-0.5, 2]));

ps_len = length(ps)
Plots.plot(ps..., size = (1200, 200 * ps_len), layout = grid(ps_len, 1), titlefont = 10)


