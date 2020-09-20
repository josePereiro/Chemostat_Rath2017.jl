
## Description
# This script produce data to be used for reporting. It uses the data generated from (./fva_pp_tINIT_models_maxent_ep.jl).

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
import Chemostat.Utils: av, va, μ, σ, bounds, norm1_stoi_err

import Chemostat_Rath2017
import Chemostat_Rath2017: RathData, DATA_KEY, load_data, save_data
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs, ecGEMs, OBJ_IDER, 
                                RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT, 
                                IDERS_MAP
const Rd = RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const ecG = ecGEMs

## Loading dat
src_file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
bundles = load_data(src_file)

## Find experimental beta
# The beta value that has the closest growth value to the experimental growth
function find_exp_β(bundle, stst)
    exp_μ = Rd.val(:μ, stst)
    βs = bundle.βs
    exp_β = βs |> first
    for β in βs
        ep_μ = av(bundle, 1, β, :ep, OBJ_IDER)
        last_ep_μ = av(bundle, 1, exp_β, :ep, OBJ_IDER)
        if abs(ep_μ - exp_μ) < abs(last_ep_μ - exp_μ)
            exp_β = β
        end
    end
    return exp_β
end

## Extracting Data
global_dict = Dict()

EXP_BETA_KEY = "exp_beta"
XIS_KEY = "xis"
BETAS_KEY = "betas"
MEAN_STOI_ERR_KEY = "mean_ep_norm_stoi_err"
MAX_STOI_ERR_KEY = "max_ep_norm_stoi_err"

##
for (model_id, dat) in bundles
    
    model_data = get!(global_dict, model_id, Dict())
    for (stst, bundle) in dat

        stst_data = get!(model_data, stst, Dict())
        exp_β = find_exp_β(bundle, stst)
        stst_data[EXP_BETA_KEY] = exp_β

        ξs = bundle.ξs
        stst_data[XIS_KEY] = ξs
        βs = bundle.βs
        stst_data[BETAS_KEY] = βs
    
        for rath_ider in RATH_IDERS_TO_PLOT
            model_ider = IDERS_MAP[rath_ider]
            for ξ in ξs
                # FBA
                key = string((ξ, :fba, rath_ider))
                stst_data[key] = (av = av(bundle, ξ, :fba, model_ider),)

                # bounds
                model = bundle[ξ, :net]
                key = string((ξ, :bounds, rath_ider))
                stst_data[key] = bounds(model, model_ider)
                
                # EP
                for β in βs
                    key = string((ξ, β, :ep, rath_ider))
                    stst_data[key] = (
                        av = av(bundle, ξ, β, :ep, model_ider),
                        va = va(bundle, ξ, β, :ep, model_ider),
                        μ = μ(bundle, ξ, β, :ep, model_ider),
                        σ = σ(bundle, ξ, β, :ep, model_ider)
                    )

                end
            end
        end # for model_ider in MODEL_IDERS_TO_PLOT


        # Stoi error
        for ξ in ξs
            metnet = bundle[ξ, :net]
            epouts = bundle[ξ, βs, :ep]
            errs = map(epouts) do epout 
                norm1_stoi_err(metnet, epout)
            end
            key = string((ξ, MEAN_STOI_ERR_KEY))
            stst_data[key] = mean.(errs)
            key = string((ξ, MAX_STOI_ERR_KEY))
            stst_data[key] = maximum.(errs)
        end

    end
end

## Save
save_data(ecG.EXTRACTED_DATA_FILE, global_dict)


