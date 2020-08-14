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

# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import Chemostat_Rath2017: DATA_KEY, Rep_Human1
const RepH1 = Rep_Human1;

using SparseArrays
using StatsBase

import Chemostat;
const Ch = Chemostat
# -

using Plots
pyplot();

using Distributions

dat_file = RepH1.COMP_FVA_HG_EXTRACTED_DATA_FILE;
dat = wload(dat_file)[DATA_KEY]
model_keys = [:orig_model, :ec_model];

## Loading Models
input_dat  = wload(RepH1.COMP_FVA_HG_INPUT_FILE);
orig_model = Ch.Utils.uncompress_model(input_dat["orig_model"]);
ec_model   = Ch.Utils.uncompress_model(input_dat["ec_model"]);

model = ec_model;
lb, ub = dat[:ec_model][:bounds];
clamp!(lb, -1000, 1000);
clamp!(ub, -1000, 1000);

function get_revs_map(model; rev_suffix = "_REV")
    M, N = size(model);
    revs_map = Dict()
    for i in 1:N
        rxn = model.rxns[i]
        if endswith(rxn, rev_suffix)
            rxn = replace(rxn, rev_suffix => "")
            revs_map[rxn] = filter(1:N) do idx
                rxn_ = model.rxns[idx]
                occursin(rxn, rxn_)
            end
        end
    end
    return revs_map
end

revs_map = get_revs_map(model);

rev_suffix = "_REV";
var_ = [];
for (i, rxn) in model.rxns |> enumerate
    if haskey(revs_map, rxn)
        idxs = revs_map[rxn]
        push!(var_, abs(ub[idxs[1]] - ub[idxs[2]]))
    elseif !endswith(rxn, rev_suffix)
        push!(var_, abs(ub[i] - lb[i]))
    end
end
var_ = filter(v -> v > zeroth, var_);

sort!(var_);
Plots.plot(var_ .+ zeroth, (1:length(var_)) ./ length(var_), xscale = :log10)

zeroth = 1e-12
infth = 3e3
p = Plots.plot(
    xlabel = "flux variability", ylabel = "Fraction of reactions", 
    xscale = :log10,
    xlim = [1e-12, 1e4]
)
bounds_sym = :bounds
for sym in model_keys
    
    lb, ub = dat[sym][bounds_sym]
    
    # cdf
    diffs_ = abs.(ub .- lb)
    nz_diffs_ = filter(diff -> zeroth < diff < infth , diffs_)
    N = length(nz_diffs_)
    xs, ys = (sort!(nz_diffs_), (1:N) / N)
    Plots.plot!(p, xs, ys, label = sym)
    
    # median
    med = median(nz_diffs_)
    Plots.scatter!(p, [med], [0.5]; 
        label = "", color = :black)
    med_str = string(round(med, digits = 6))
    Plots.annotate!(p, [(med, 0.5, text(med_str, 10, :right, :top))])
end
p







diff_cdf_sym = :diff_pdf
p = Plots.plot(
    xlabel = "flux variability", ylabel = "Fraction of reactions", 
    xscale = :log10,
    xlim = [1e-12, 1e4]
)
for sym in model_keys
    xs, ys, mean, std = dat[sym][diff_cdf_sym];
    @show mean
    Plots.plot!(p, xs, ys, label = sym)
#     Plots.scatter!(p, [mean], [0.5], xerr = [std], label = "", color = :black)
end
Plots.scatter!(p, [], [], label = "mean/std", color = :black)
p




