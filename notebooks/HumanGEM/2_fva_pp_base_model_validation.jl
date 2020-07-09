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

# +
using DataFrames
using Serialization
using Dates
using StatsBase

using Plots
pyplot()
# -

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
H1 = Chemostat_Rath2017.Human1

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# ## Model

model = deserialize(H1.FVA_PP_BASE_MODEL_FILE);
obj_ider = "biomass_human"
obj_idx = Ch.Utils.rxnindex(model, obj_ider);

H1.load_base_intake_info();

# +
model = deserialize(H1.FVA_PP_BASE_MODEL_FILE)
# foreach(Ch.Utils.exchanges(model)) do exch
#     # "HMR_9034" ex_glc
#     Ch.Utils.lb!(model, exch, -1000)
# end

H1.load_base_intake_info()
for stst in Rd.ststs |> sort
    ξ = Rd.val(:ξ, stst)
    intake_info = H1.stst_base_intake_info(stst)
    
    Ch.SteadyState.apply_bound!(model, ξ, intake_info);
    fba_av = Ch.LP.fba(model, obj_ider).obj_val
    exp_av = Rd.val(:μ, stst)
    println(stst)
    println("xi: $ξ")
    println("growth: fba: $fba_av, exp: $exp_av")
end
# -

# ### FBA

ξs = range(10, 1e5, length = 25)
ξs = [ξs; Rd.val("ξ", Rd.ststs)] |> collect |> unique |> sort;

# +
model = deserialize(H1.FVA_PP_BASE_MODEL_FILE);
obj_ider = "biomass_human"
obj_idx = Ch.Utils.rxnindex(model, obj_ider);

fbaouts = Dict()

for (ststi, stst) in Rd.ststs |> enumerate
    
    # intake info
    intake_info = Dict()
    delete!(intake_info, "HMR_9136") # Not present in the model
    delete!(model.intake_info, "HMR_9136")

    fbaouts[stst] = []
    for (ξi, ξ) in ξs |> enumerate
        
        print("Doing stst[$ststi/$(length(Rd.ststs))]: $stst, xi[$ξi/$(length(ξs))]: $ξ                    \r")
        flush(stdout)
        
        Ch.SteadyState.apply_bound!(model, ξ, intake_info)
        fbaout = Ch.LP.fba(model, obj_ider);
#         fbaout = fba(model, obj_ider);
        
        push!(fbaouts[stst], fbaout)
    end
end
println("Done!!!"," "^150)
# -

# ## Plots

colors = Dict()
_colors = Rd.ststs |> length |> Plots.distinguishable_colors
for (stst, color) in zip(Rd.ststs, _colors)
    colors[stst] = color
end

p = plot()
for (stst, fbas) in fbaouts
    avs = [Ch.Utils.av(model, fba, obj_ider) for fba in fbas]
    plot!(ξs, avs, color = colors[stst], label = "")
    
    expξ = Rd.val(:ξ, stst)
    expμ = Rd.val(:μ, stst)
    scatter!([expξ], [expμ], color = colors[stst], label = stst, m = 10)
end
p


