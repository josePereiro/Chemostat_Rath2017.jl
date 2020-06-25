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
import DataFrames: DataFrame
import MAT
import CSV
import Serialization: serialize, deserialize

using Plots
pyplot();

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
M = Chemostat_Rath2017.MODEL1105100000
Rd = Chemostat_Rath2017.RathData
# -

notebook_name = "medium_profile"

# +
# GLC
stst = "A"
ξ = Rd.val("ξ", stst)
μs = Dict()
base_model = deserialize(M.BASE_MODEL_FILE);
obj_ider = "BIOMASS";

conc_fs = 0.0:0.1:1.0
intake_info = M.stst_base_intake_info(stst)

for (intake, info) in intake_info
    
    conc0 = info["c"]
    μs[intake] = []
    
    for (i, conc_f) in conc_fs |> enumerate

        conc = conc0 * conc_f
        intake_info[intake]["c"] = conc


        Ch.SteadyState.apply_bound!(base_model, ξ, intake_info)
        
        μ = 0.0
        try
            fbaout = Ch.FBA.fba(base_model, obj_ider)
            μ = Ch.Utils.av(base_model, fbaout, obj_ider)
        catch err
            err isa InterruptException && rethrow(err)
        end
        push!(μs[intake], μ)
        
        # Progress
        print("Done $intake conc[$i/ $(length(conc_fs))]: $conc mM , μ: $μ            \r"); flush(stdout)
        
    end
    
    intake_info[intake]["c"] = conc0
end
println("Done!!!                                                        ")

# +
## Relevants
glc_intake = M.exch_met_map[M.mets_map["GLC"]]
th = 0.2 # Selection threshold
rel_μs = filter(μs) do p
    intake, μs_ = p
    μs_ = μs_[2:end] # exclude conc = 0
    return maximum(μs_) - minimum(μs_) >= th * maximum(μs_) || intake == glc_intake
end

rel_ps = []
for (intake, μs_) in rel_μs
    met = M.exch_met_map[intake]
    p = plot(title = met, xticks = false, yticks = false, ylim = [0.0, maximum(μs_)*1.1])
    plot!(conc_fs, μs_, label = "")
    push!(rel_ps, p)
end
# -

n = length(rel_ps)
c = 4
r = floor(Int, n/c)
r = (n % c) == 0 ? r : r += 1
println("grid: ", (r, c))
plot(rel_ps..., titlefont = 10, size = [400, 110])

filename = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(notebook_name).png");
savefig(p, filename)
println(relpath(filename), " created!!!")
