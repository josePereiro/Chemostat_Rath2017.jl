## ---------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")
import Weave: weave
import Dates: now

## ---------------------------------------------------------------
jmd_dat = []

## ---------------------------------------------------------------
# Weaveing: MARKDOWN
push!(jmd_dat, 
"""
---
title: FIR filter design with Julia
author: Jose
date: $(now())
---
"""
)

## ---------------------------------------------------------------
# Import
push!(jmd_dat,  
# ```julia; echo = false; results = "hidden"
"""
```julia
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

import SparseArrays
import Distributions: mean
import Chemostat_Rath2017: Chemostat
import Chemostat_Rath2017.Chemostat.LP: MathProgBase
import Chemostat_Rath2017.Human1: HumanGEM, RathData, ecGEMs
const ChU = Chemostat.Utils
const HG = HumanGEM
const Rd = RathData
const ecG = ecGEMs

using Plots

maxent_dat = ChU.load_data(ecG.MAXENT_FBA_EB_BOUNDLES_FILE)
commit_hash = ChU.load_commit_short_hash(ecG.MAXENT_FBA_EB_BOUNDLES_FILE)
fig_title = string(nameof(ecG), " [", commit_hash, "]") 

color_pool = [:orange, :blue, :red, :black, :violet, 
:gray, :green, :brown, :magenta];
colors = Dict(exp => rand(color_pool) for exp in Rd.exps)

closest_βs = Dict()
for (model_id, bundles) in maxent_dat
    closest_βs[model_id] = Dict()
    for (exp, bundle) in bundles
        exp_xi = Rd.val("ξ", exp)
        exp_mu = Rd.val("μ", exp)
        closest_βs[model_id][exp] = ChU.find_closest_beta(bundle, exp_xi, exp_mu, HG.OBJ_IDER)
    end
end

function growth_vs_beta()
    eps_ = 1e-5
    p = plot(title = fig_title, 
        xlabel = "log beta", ylabel = "log growth rate", 
        legend = false
    )
    fun1 = log10 #(x) -> x
    for (model_id, bundles) in maxent_dat
        sbundles = sort(collect(bundles); by = (x) -> x[1])
        for (exp, bundle) in sbundles
            # exp != 7 && continue
            exp_growth = Rd.val("μ", exp)
            plot!(p, fun1.(bundle.βs), fill(fun1(exp_growth), length(bundle.βs)); 
                ls = :dash, color = colors[exp], lw = 3, label = "")
            
            exp_xi = Rd.val("ξ", exp)
            fba_growth = ChU.av(bundle, exp_xi, :fba, HG.OBJ_IDER)
            fba_growth += eps_
            plot!(p, fun1.(bundle.βs), fill(fun1(fba_growth), length(bundle.βs)); 
                ls = :dot, color = colors[exp], lw = 3, label = "")

            ep_growths = ChU.av(bundle, exp_xi, bundle.βs, :ep, HG.OBJ_IDER)
            ep_growths .+= eps_
            ep_stds = sqrt.(ChU.va(bundle, exp_xi, bundle.βs, :ep, HG.OBJ_IDER))
            plot!(p, fun1.(bundle.βs), fun1.(ep_growths); 
                label = exp, lw = 3, color = colors[exp])

        end
    end
    return p
end

growth_vs_beta()
```
"""
)


## ---------------------------------------------------------------
# weaveing
jmd_file = joinpath(@__DIR__, "make_report.jmd")
write(jmd_file, join(jmd_dat, "\n"))
weave(jmd_file)