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

using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using DataFrames
using Serialization
using StatsBase
using Statistics
using Measures
using Plots
pyplot()

import Chemostat
Ch = Chemostat;

using Revise
import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
M = Chemostat_Rath2017.MODEL1105100000;

# ---
# ## Meta
# ---

notebook_name = "fva_pp_base_model_plot";
data_notebook_name = "fva_pp_base_model_maxent_ep_epsconv_study";

# ---
# ## MaxEnt-EP
# --- 

# Check cached objects
cache_file = joinpath(M.MODEL_CACHE_DATA_DIR, "$(data_notebook_name)___boundles.jls");
params, data = deserialize(cache_file)
println("loaded $(relpath(cache_file))!!!")

boundles = Dict()
for (state, boundle) in data
    boundles[state] = boundle
end

# +
# ### Cleaning boundles, keep only the complete errorless indexes
# function get_clean(boundle)
#     βs_ = sort(boundle.βs, rev = true)
#     ξs_ = sort(boundle.ξs)
    
#     errorless_βi = trues(length(βs_))
#     for ξ in ξs_
#         for (βi, β) in enumerate(βs_)
#             if haskey(boundle, ξ, β, :ep)
#                 data = Ch.Utils.get_data(boundle, ξ, β, :ep)
#                 errorless_βi[βi] = data isa Ch.Utils.EPout && errorless_βi[βi]
#             else
#                 errorless_βi[βi] = false
#             end
#         end
#     end
#     all(errorless_βi .== false) && error("Boundle is too dirty!!!")
    
#     errorless_βi = βs_[errorless_βi]
#     errorless_boundle = Ch.Utils.ChstatBoundle()
#     for ξ in ξs_
#         model = Ch.Utils.get_data(boundle, ξ, :net)
#         fbaout = Ch.Utils.get_data(boundle, ξ, :fba)
#         Ch.Utils.add_data!(errorless_boundle, ξ, :net, model)
#         Ch.Utils.add_data!(errorless_boundle, ξ, :fba, fbaout)
        
#         for β in errorless_βi
#             epout = Ch.Utils.get_data(boundle, ξ, β, :ep)
#             Ch.Utils.add_data!(errorless_boundle, ξ, β, :ep, epout)
#         end
#     end
#     return errorless_boundle
# end

# +
# for (stst, boundle) in boundles
#     boundles[stst] = get_clean(boundle)
# end
# # This three ststs are equivalents (repetitions of the same experiment)
# # boundles["A"] = boundles["C"]
# # boundles["B"] = boundles["C"];
# -

obj_ider = params["obj_ider"]; # all models are equals in this sense
ep_epsconvs = params["ep_epsconvs"];

# ---
# ## Plots
# ---

function commom_data(stst, ep_epsconv, model_ider, idxs...)
    
    boundle = boundles[(stst, ep_epsconv)]
    
    exp_ider = model_ider == obj_ider ? "μ" : "q" * M.mets_map[M.exch_met_map[model_ider]]
    sense = model_ider == obj_ider ? 1 : -1

    exp_av = sense * Rd.val(exp_ider, stst)
    exp_av_err = Rd.err(exp_ider, stst)
    
    ep_av = Ch.Utils.av(boundle, idxs..., model_ider)
    ep_std = sqrt(Ch.Utils.va(boundle, idxs..., model_ider))  
    
    return exp_av, exp_av_err, ep_av, ep_std
end

exchs = [M.exch_met_map[M.mets_map[met]] for met in Rd.msd_mets];

stst_colors = Dict()
for (stst, c) in zip(Rd.ststs, Plots.distinguishable_colors(length(Rd.ststs)))
    stst_colors[stst] = c
end
stst_colors["A"] = stst_colors["C"];
stst_colors["B"] = stst_colors["C"];

epsconv_colors = Dict()
for (eps, c) in zip(ep_epsconvs, Plots.distinguishable_colors(length(ep_epsconvs)))
    epsconv_colors[eps] = c
end

Plots.plot(ps_iders..., layout = grid(1, length(ps_iders)), 
    size = [4000,2500], titlefont = 10, margin = 5mm)

# ### Stoi err

# +
p = Plots.plot()
stst_ = Rd.ststs
for stst in ststs_
    
end
