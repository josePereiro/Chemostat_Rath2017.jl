using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using Dates
using StatsBase
using DataFrames
using SparseArrays
using MathProgBase
using Serialization
using Distributions

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.Utils: av, va, μ, σ, bounds, is_exchange,
                        norm1_stoi_err, norm2_stoi_err,
                        find_closest_beta, metindex, met_rxns

import UtilsJL: save_data, load_data, to_symbol_dict
import Chemostat_Rath2017
import Chemostat_Rath2017: DATA_KEY, RathData
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs, ecGEMs, OBJ_IDER, 
                                RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT, 
                                IDERS_MAP, RATH_GROWTH_IDER
const Rd = RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const ecG = ecGEMs

## Loading data
src_file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
bundles = load_data(src_file)

## Dev data
model_id = "GTEx-brain"
stst = "E"
bundle = bundles[model_id][stst]
exp_ξ = Rd.val(:ξ, stst)
exp_μ = Rd.val(:D, stst)
model = bundle[exp_ξ, :net]
fbaout = bundle[exp_ξ, :fba]
exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, OBJ_IDER)
epout = bundle[exp_ξ, exp_β, :ep];


##
using StatsPlots
# gr(size = (600, 500))

##

## Total Correlations
function total_correlation()
    lim = 0.2
    lims = [-lim, lim]
    fbap = plot(
        ylim = lims, 
        xlim = lims, 
        title = "FBA correlation", 
        xlabel = "exp", ylabel = "model"
    )
    lim = 5
    lims = [-1, lim]
    epp = plot(
        ylim = lims, 
        # xlim = lims, 
        title = "EP correlation", 
        xlabel = "exp", ylabel = "model"
    )
    for (model_id, stst_dat) in bundles
        for (stst, bundle) in stst_dat
            exp_ξ =  Rd.val(:ξ, stst)
            exp_μ =  Rd.val(:D, stst)
            exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, OBJ_IDER)

            for rider in RATH_IDERS_TO_PLOT
                sense = rider == RATH_GROWTH_IDER ? 1 : -1 # TODO: package this
                mider = IDERS_MAP[rider]
                exp_av = sense * Rd.qval(rider, stst)
                exp_err = Rd.qerr(rider, stst)

                # FBA
                mod_av = av(bundle, exp_ξ, :fba, mider)
                scatter!(fbap, [exp_av], [mod_av]; xerr = [exp_err], label = "")

                # EP
                mod_av = av(bundle, exp_ξ, exp_β,:ep, mider)
                mod_err = sqrt(va(bundle, exp_ξ, exp_β,:ep, mider))
                scatter!(epp, [exp_av], [mod_av]; xerr = [exp_err], yerr = [mod_err], label = "")
            end
        end
    end
    return plot([fbap, epp]..., size = [800, 400])
end
total_correlation()

## stoi_err_by_compartments
# errs = log.(norm2_stoi_err(model, epout) .+ 1e-3);
# histogram(errs; normalize = true, xlabel = "log err", ylabel = "prob dens")

function stoi_err_by_compartments(get_comp::Function, model, out)
    
    errs = norm2_stoi_err(model, out)

    ## Split error by compartment
    cerrs = Dict()
    for (meti, met) in model.mets |> enumerate
        c = get_comp(model, met) 
        push!(get!(cerrs, c, []), errs[meti])
    end

    ps = []
    xlims_ = [-6,4]
    bins_ = 100
    for (c, errs) in cerrs
        log_err = log10.(errs .+ 1e-5) 
        p = histogram(log_err; label = string(c), 
            xlabel = "log err", ylabel = "prob", 
            xlims = xlims_, bins = bins_,
            normalize = true)
        push!(ps, p)
    end
    ##
    plot(ps..., size = [1000, 400])

end

cerrs = stoi_err_by_compartments(model, epout) do model, met
    rxns = met_rxns(model, met)
    return any(is_exchange(model, rxn) for rxn in rxns) ? "ext" : "int"
end

##
errs = norm1_stoi_err(model, epout);
maximum

##
plot(["A", "B"], [2,3])