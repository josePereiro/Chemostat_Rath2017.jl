# # -*- coding: utf-8 -*-
# # ---
# # jupyter:
# #   jupytext:
# #     cell_metadata_filter: -all
# #     formats: jl
# #     text_representation:
# #       extension: .jl
# #       format_name: light
# #       format_version: '1.5'
# #       jupytext_version: 1.3.2
# #   kernelspec:
# #     display_name: Julia 1.1.0
# #     language: julia
# #     name: julia-1.1
# # ---

# # +
# using DataFrames
# using Serialization
# using Dates
# using StatsBase

# using Plots
# pyplot()
# # -

# # run add "https://github.com/josePereiro/Chemostat" in the 
# # julia Pkg REPL for installing the package
# import Chemostat
# Ch = Chemostat

# import Chemostat_Rath2017
# Rd = Chemostat_Rath2017.RathData
# HG = Chemostat_Rath2017.HumanGEM

# # This just check that the script is run in the
# # package enviroment
# Chemostat_Rath2017.check_env();

# # ## Cached files

# # Chached files
# files = ["fva_pp_base_model_maxent_ep_v1___temp_cache___state_12262029076231720947.jls", 
#     "fva_pp_base_model_maxent_ep_v1___temp_cache___state_13145855622767250497.jls"]
# files = joinpath.("/Volumes/Store/Temp Results/Chemostat_Rath2017/cache", files)
# @assert isfile.(files) |> all
# data = deserialize.(files);

# bundles = Dict()
# for dat in data
#     stst, bundle = dat
#     bundles[stst] = bundle
# end

# obj_ider = "biomass_human"
# good_βs = Dict()
# for (stst, bundle) in bundles
#     good_βs[stst] = filter((β) -> bundle[1, β, :ep] isa Ch.Utils.EPout, bundle.βs);
# end

# # +
# colors = [:blue, :red]#Plots.distinguishable_colors(length(bundles))
# p = plot(title = "HumanGEM", xlabel = "beta", ylabel = "alpha")
# for (i, (stst, bundle)) in bundles |> enumerate
#     ep_avs = Ch.Utils.av(bundle, 1, good_βs[stst], :ep, obj_ider);
#     fba_av = Ch.Utils.av(bundle, 1, :fba, obj_ider)
#     scatter!(good_βs[stst], ep_avs, color = colors[i], label = "", m = 10)
#     hline!([fba_av], color = colors[i], lw = 3, label = "")
#     hline!([Rd.val(:μ, stst)], color = colors[i], lw = 3, label = "", ls = :dash)
#     vline!([maximum(bundle.βs)], label = "", color = colors[i])
#     scatter!([],[], label = stst, marker = :square, ms = 30, color = colors[i])
# end
# scatter!([],[], ls = :dash, color = :black, label = "ep")
# plot!([],[], ls = :dash, color = :black, label = "exp")
# plot!([],[], color = :black, label = "fba")

# p
# # -

# ## Stoi err
# p = plot(title = "HumanGEM", xlabel = "beta", ylabel = "abs norm stoi error")
# for (i, (stst, bundle)) in bundles |> enumerate
#     βs = good_βs[stst] |> sort
#     metnet = bundle[1, :net]
#     epouts = bundle[1, βs, :ep]
#     errs = map(epouts) do epout 
#         Ch.Utils.norm1_stoi_err(metnet, epout)
#     end
#     plot!(βs, mean.(errs), label = "", color = colors[i], lw = 2, ls = :dash)
#     plot!(βs, maximum.(errs), label = "", color = colors[i], lw = 2)
#     scatter!([],[], label = stst, marker = :square, ms = 30, color = colors[i])
# end
# plot!([],[], ls = :dash, color = :black, label = "mean")
# plot!([],[], color = :black, label = "max")
# p


