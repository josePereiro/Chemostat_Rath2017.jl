import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## --------------------------------------------------------------------
using DataFrames
using Serialization
using Dates
using StatsBase
using JSON
using SparseArrays
using MathProgBase

## --------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.Utils: av, va, μ, σ, bounds, ChstatBundle, rxn_mets,
                        norm2_stoi_err, stoi_err, met_rxns, lb, ub, rxn_str,
                        summary, is_exchange, similar_rxns, exchanges, 
                        uncompressed_copy, MetNet, rxnindex
import Chemostat.SteadyState: apply_bound!
import Chemostat.LP: fva_preprocess, fva

import Chemostat_Rath2017
import Chemostat_Rath2017: RathData, DATA_KEY, load_data, save_data
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs, ecGEMs, OBJ_IDER, 
                                RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT, 
                                IDERS_MAP, RATH_GROWTH_IDER
const Rd = RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const ecG = ecGEMs

## --------------------------------------------------------------------
using Plots

##
function find_closest_β(bundle, target, ider)
    βs = bundle.βs
    exp_β = βs |> first
    for β in βs
        ep_μ = av(bundle, 1, β, :ep, ider)
        last_ep_μ = av(bundle, 1, exp_β, :ep, ider)
        if abs(ep_μ - target) < abs(last_ep_μ - target)
            exp_β = β
        end
    end
    return exp_β
end

## Loading dat
src_file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
bundles = load_data(src_file);

##
model_id = "GTEx-brain"
stst = "E"
bundle = bundles[model_id][stst];

##
exp_μ = Rd.val(:μ, stst)
exp_β = find_closest_β(bundle, exp_μ, OBJ_IDER)
xi = 1
beta = exp_β
model = bundle[xi, :net]
fbaout = bundle[xi, :fba]
epout = bundle[xi, beta, :ep];

##
## COMPARE FBA -- EP
rf(n) = round(n; digits = 3)
intake_info = HG.stst_base_intake_info(stst)
th = 50
exchs = exchanges(model)
sep = ","
csv = [join(["id", "fbaflx", "epflx", "epflx_va", "lb", "ub", "eq", "imp"], sep)]
for rxni in exchs
    rxn = model.rxns[rxni]
    important = rxn in MODEL_IDERS_TO_PLOT || haskey(intake_info, rxn)
    rxn_str_ = important ? rxn * "*" : rxn
    epflx = av(bundle, xi, beta, :ep, rxn)
    epflx_va = va(bundle, xi, beta, :ep, rxn)
    fbaflx = rf(av(bundle, xi, :fba, rxn))
    lb_ = lb(model, rxn)
    ub_ = ub(model, rxn)
    eq = rxn_str(model, rxn)
    push!(csv, join(string.([rxn_str_, fbaflx, epflx, epflx_va, lb_, ub_, eq, important ? 1 : 0]), sep))
    !important && !(abs(epflx) > th || abs(fbaflx) > th) && continue
    println(rpad(string(rxn_str_, " "), 30, "-"), " ", 
        rpad(string("b: ", (rf(lb_), rf(ub_)), " "), 25, "-"), " ",
        rpad(string("fba: ", rf(fbaflx), " "), 15, "-"), " ",
        rpad(string("ep: ", rf(epflx), " "), 15, "-"), " ",
        rpad(string("ep_va: ", rf(epflx_va)), 10)
    )
end
# 0.298
open("dat.csv", "w") do io
    for line in csv
        println(io, line)
    end
end

##
eps = 1e-5
epflx = abs.(av(bundle, xi, beta, :ep, exchs)) .+ eps 
epflx_va = abs.(va(bundle, xi, beta, :ep, exchs)) .+ eps;
scatter(epflx, epflx_va; xscale = :log10, yscale = :log10, label = "", 
    xlabel = "flx", ylabel = "var")

##
# rxn = "HMR_9133" 
met_ = "m02819s"
for rxni in met_rxns(model, met_)
    rxn = model.rxns[rxni]
    summary(model, rxn)
end

##
similar_dict = Dict()
for rxn in MODEL_IDERS_TO_PLOT
    !is_exchange(model, rxn) && continue
    met = rxn_mets(model, rxn) |> first
    similars = vcat(similar_rxns(model, met_rxns(model, met); verbose = false)...);
    similar_dict[rxn] = length(similars)
end


##
iders = similar_dict |> keys |> collect
vals = similar_dict |> values |> collect
bar(iders, vals; size = [800, 400])

##
range_dict = Dict()
for rxn in MODEL_IDERS_TO_PLOT
    range_dict[rxn] = abs(lb(model, rxn) - ub(model, rxn))
end
iders = range_dict |> keys |> collect
vals = range_dict |> values |> collect
bar(iders, vals; size = [800, 400])

##
intake_info = HG.stst_base_intake_info(stst)
apply_bound!(model, xi, intake_info; 
        emptyfirst = true, ignore_miss = true)

##
for (rxn, dat) in model.intake_info
    metis = rxn_mets(model, rxn)
    mets = model.mets[metis]
    names = [get!(HG.readable_met_ids_map, met, "") for met in mets]
    brange = abs(lb(model, rxn) - ub(model, rxn))
    # brange > 100 && continue
    println(rxn, " [", names, "]: ", lb(model, rxn))
    
end

##
similars = similar_rxns(model; verbose = false);

##
# Exchanges
exchs = exchanges(model)
for exchi in exchs
    exch = model.rxns[exchi]
    flx = av(model, epout, exch)
    flx >= 0.0 && continue
    println(exch, " flx:", flx, " ----- ", rxn_str(model, exch))
end


## FVA PREPROCESSING
fva_temp_file = "fva_res.bson"
exchs = exchanges(model)
model_dict = struct2dict(model) |> uncompressed_copy
fvalb, fvaub = fva(MetNet(;model_dict...), 
    exchs; verbose = true);
save_data(fva_temp_file, (fvalb, fvaub))

##
fvalb, fvaub = load_data(fva_temp_file)
##
println()
th = 100
for (i, rxni) in exchs |> enumerate
    rxn = model.rxns[rxni]
    lb, ub = fvalb[i], fvaub[i]
    !(abs(lb) > th || abs(ub) > th) && continue
    println(rpad(rxn, 30), " ", 
        rpad(string("lb: ", lb), 35),
        rpad(string("ub: ", ub), 25)
    )
end
