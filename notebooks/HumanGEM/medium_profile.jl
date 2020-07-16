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

using Revise

import MathProgBase.HighLevelInterface: linprog
import Clp: ClpSolver


betas = 10.0 .^ range(2, 4, length = 20) |> collect |> reverse
betas[15]

Matrix{Float64}(rand(10,10))



30778.810222077907 |> log10

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
HG = Chemostat_Rath2017.HumanGEM
Rd = Chemostat_Rath2017.RathData
# -

notebook_name = "medium_profile";

Ch.LP.

all(uvals .== toy_model.ub) && all(lvals .== toy_model.lb)





# +
# GLC
stst = "A"
ξ = Rd.val("ξ", stst)
μs = Dict()
base_model = deserialize(HG.BASE_MODEL_FILE);
obj_ider = "biomass_human";

conc_fs = 0.0:0.01:1.0
intake_info = HG.stst_base_intake_info(stst)
errs = []

for (intake, info) in intake_info
    
    conc0 = info["c"]
    μs[intake] = []
    
    for (i, conc_f) in conc_fs |> enumerate

        conc = conc0 * conc_f
        intake_info[intake]["c"] = conc


        Ch.SteadyState.apply_bound!(base_model, ξ, intake_info)
        
        μ = 0.0
        try
            fbaout = Ch.LP.fba(base_model, obj_ider)
            μ = Ch.Utils.av(base_model, fbaout, obj_ider)
        catch err
            err isa InterruptException && rethrow(err)
            push!(errs, "Error at intake: $intake, conc_f: $conc_f")
        end
        push!(μs[intake], μ)
        
        # Progress
        print("Done $intake conc[$i/ $(length(conc_fs))]: $conc mM , μ: $μ            \r"); flush(stdout)
        
    end
    
    intake_info[intake]["c"] = conc0
end
println("Done!!!                                                        ")
println("Errors")
println.(errs)

# +
## Relevants
glc_intake = HG.exch_met_map[HG.mets_map["GLC"]]
th = 0.2 # Selection threshold
rel_μs = filter(μs) do p
    intake, μs_ = p
    μs_ = μs_[2:end] # exclude conc = 0
    return maximum(μs_) - minimum(μs_) >= th * maximum(μs_) || intake == glc_intake
end
rel_intakes = rel_μs |> keys |> collect

rel_ps = []
for (intake, μs_) in rel_μs
    readable_met = HG.readable_met_ids_map[HG.exch_met_map[intake]]
    p = plot(title = readable_met, xticks = false, yticks = false, ylim = [0.0, maximum(μs_)*1.1])
    plot!(conc_fs, μs_, label = "")
    push!(rel_ps, p)
end
file_ = joinpath(HG.MODEL_PROCESSED_DATA_DIR, "$(notebook_name).jls")
serialize(file_, (conc_fs, μs)) # Cache
println(relpath(file_), " created!!!")
# -

n = length(rel_ps)
c = 4
r = floor(Int, n/c)
r = (n % c) == 0 ? r : r += 1
println("grid: ", (r, c))
empty_plot =  plot(legend=false, grid=false, foreground_color_subplot=:white)
while length(rel_ps) < r*c
    push!(rel_ps, empty_plot)
end
p = plot(rel_ps..., layout = grid(r, c), size = [r * 100, c * 150], titlefont = 10)

filename = joinpath(HG.MODEL_FIGURES_DATA_DIR, "$(notebook_name).png");
savefig(p, filename)
println(relpath(filename), " created!!!")


