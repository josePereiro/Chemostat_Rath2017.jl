
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
import Chemostat.Utils: av, va, μ, σ, bounds, norm_abs_stoi_err

import Chemostat_Rath2017
import Chemostat_Rath2017: RathData, DATA_KEY
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs, ecGEMs, OBJ_IDER, 
                                RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT, 
                                IDERS_MAP
const Rd = RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const ecG = ecGEMs

## Loading dat
src_file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
boundles = wload(src_file)[DATA_KEY];
println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")

## Find experimental beta
# The beta value that has the closest growth value to the experimental growth
function find_exp_β(boundle, stst)
    exp_μ = Rd.val(:μ, stst)
    βs = boundle.βs
    exp_β = βs |> first
    for β in βs
        ep_μ = av(boundle, 1, β, :ep, OBJ_IDER)
        last_ep_μ = av(boundle, 1, exp_β, :ep, OBJ_IDER)
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
for (model_id, dat) in boundles
    
    model_data = get!(global_dict, model_id, Dict())
    for (stst, boundle) in dat

        stst_data = get!(model_data, stst, Dict())
        exp_β = find_exp_β(boundle, stst)
        stst_data[EXP_BETA_KEY] = exp_β

        ξs = boundle.ξs
        stst_data[XIS_KEY] = ξs
        βs = boundle.βs
        stst_data[BETAS_KEY] = βs
    
        for rath_ider in RATH_IDERS_TO_PLOT
            model_ider = IDERS_MAP[rath_ider]
            for ξ in ξs
                # FBA
                key = string((ξ, :fba, rath_ider))
                stst_data[key] = (av = av(boundle, ξ, :fba, model_ider),)

                # bounds
                model = boundle[ξ, :net]
                key = string((ξ, :bounds, rath_ider))
                stst_data[key] = bounds(model, model_ider)
                
                # EP
                for β in βs
                    key = string((ξ, β, :ep, rath_ider))
                    stst_data[key] = (
                        av = av(boundle, ξ, β, :ep, model_ider),
                        va = va(boundle, ξ, β, :ep, model_ider),
                        μ = μ(boundle, ξ, β, :ep, model_ider),
                        σ = σ(boundle, ξ, β, :ep, model_ider)
                    )

                end
            end
        end # for model_ider in MODEL_IDERS_TO_PLOT


        # Stoi error
        for ξ in ξs
            metnet = boundle[ξ, :net]
            epouts = boundle[ξ, βs, :ep]
            errs = map(epouts) do epout 
                norm_abs_stoi_err(metnet, epout)
            end
            key = string((ξ, MEAN_STOI_ERR_KEY))
            stst_data[key] = mean.(errs)
            key = string((ξ, MAX_STOI_ERR_KEY))
            stst_data[key] = maximum.(errs)
        end

    end
end

## Save
save_file = ecG.EXTRACTED_DATA_FILE
tagsave(save_file, Dict(DATA_KEY => global_dict))
println(relpath(save_file), " created!!!, size: ", filesize(save_file), " bytes")

