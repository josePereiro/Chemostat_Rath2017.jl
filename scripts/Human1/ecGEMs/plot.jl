using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using DataFrames
using Serialization
using Dates
using StatsBase
using JSON
using SparseArrays
using MathProgBase

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.Utils: av, va, μ, σ, bounds, norm_abs_stoi_err

import Chemostat_Rath2017
import Chemostat_Rath2017: RathData, DATA_KEY
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs, ecGEMs, OBJ_IDER, 
                                RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT, 
                                IDERS_MAP, RATH_GROWTH_IDER
const Rd = RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const ecG = ecGEMs

## Loading data
src_file = ecG.EXTRACTED_DATA_FILE
extracted_dat = wload(src_file)[DATA_KEY]
println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")

EXP_BETA_KEY = "exp_beta"
XIS_KEY = "xis"
BETAS_KEY = "betas"
MEAN_STOI_ERR_KEY = "mean_ep_norm_stoi_err"
MAX_STOI_ERR_KEY = "max_ep_norm_stoi_err"

## Plots
using Plots
pyplot();

## Total Correlations
lim = 0.2
lims = [-lim, lim]
fbap = plot(
    ylim = lims, 
    xlim = lims, 
    title = "FBA correlation", 
    xlabel = "exp", ylabel = "model"
)
lim = 1
lims = [-lim, lim]
epp = plot(
    ylim = lims, 
    # xlim = lims, 
    title = "EP correlation", 
    xlabel = "exp", ylabel = "model"
)
for (model_id, model_dat) in extracted_dat
    for (stst, dat) in model_dat
        exp_ξ =  Rd.val(:ξ, stst)
        exp_β = dat[EXP_BETA_KEY]

        for rider in RATH_IDERS_TO_PLOT
            sense = rider == RATH_GROWTH_IDER ? 1 : -1 # TODO: package this
            mider = IDERS_MAP[rider]
            exp_av = sense * Rd.qval(rider, stst)
            exp_err = Rd.qerr(rider, stst)

            # FBA
            key = string((exp_ξ, :fba, rider))
            mod_av = dat[key].av
            scatter!(fbap, [exp_av], [mod_av]; xerr = [exp_err], label = "")

            # EP
            key = string((exp_ξ, exp_β, :ep, rider))
            mod_av = dat[key].av
            mod_err = sqrt(dat[key].va)
            scatter!(epp, [exp_av], [mod_av]; xerr = [exp_err], yerr = [mod_err], label = "")
        end
    end
end
plot([fbap, epp]..., size = [800, 400])

