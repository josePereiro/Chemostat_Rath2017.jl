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

# ---
# ## Description

# This script produce data to be used for reporting. It uses the data generated from [fva_pp_tINIT_models_maxent_ep](./fva_pp_tINIT_models_maxent_ep.jl).

using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using DataFrames
using Serialization
using Dates
using StatsBase
using JSON

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

import Chemostat_Rath2017
const ChR = Chemostat_Rath2017
const Rd = ChR.RathData
const HG = ChR.HumanGEM
const tIG = ChR.tINIT_GEMs;

# ## Tools

function load_boundles(dat_file)
    dat = deserialize(dat_file)
    dat_id = dat.id
    boundles = Dict()
    for dat in dat.res
        stst, boundle = dat
        boundles[stst] = boundle
    end
    return boundles
end

# The current maxent_ep algorithm may return error data points
function get_good_βs(boundles)
    good_βs = Dict()
    for (stst, boundle) in boundles
        good_βs[stst] = filter(boundle.βs) do β
            dat = boundle[1, β, :ep]
            dat isa Ch.Utils.EPout && dat.status == :converged
        end
    end
    return good_βs
end
# println("good betas: ", length(good_βs))

# +
function find_exp_β(boundles, good_βs)
    exp_βs = Dict()
    for (stst, boundle) in boundles
        exp_μ = Rd.val(:μ, stst)
        βs = good_βs[stst]
        exp_β = βs |> first
        for β in βs
            ep_μ = Ch.Utils.av(boundle, 1, β, :ep, obj_ider)
            last_ep_μ = Ch.Utils.av(boundle, 1, exp_β, :ep, obj_ider)
            if abs(ep_μ - exp_μ) < abs(last_ep_μ - exp_μ)
                exp_β = β
            end
        end
        exp_βs[stst] = exp_β
    end
    return exp_βs
end

# println("\nExp_beta")
# ep_μ = Ch.Utils.av(boundle, 1, exp_β, :ep, obj_ider)
# println("stst: ", stst)
# println("good_betas: ", length(βs))
# println("exp beta: ", exp_β)
# println("exp μ: ", exp_μ)
# println("ep μ:  ", ep_μ)
# println()
# -

# ## MaxEnt EP results

dat_files = filter(readdir(tIG.MODEL_PROCESSED_DATA_DIR)) do file
    startswith(file, "fva_pp_tINIT_models_maxent_ep___") && endswith(file, "___boundles.jls")
end
println("Found: ")
dat_files .|> println
dat_files = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, dat_files);
obj_ider = "biomass_human";

# ## Extract Data

for dat_file in dat_files
    
    boundles = load_boundles(dat_file)
    good_βs = get_good_βs(boundles)
    exp_β = find_exp_β(boundles, good_βs)
    data = Dict()
    
    for (stst, boundle) in boundles
        
        data[stst] = Dict()
        data[stst]["βs"] = good_βs[stst]
        data[stst]["ξs"] = boundle.ξs
        data[stst]["exp_β"] = exp_β[stst]
        data[stst]["ider_map"] = Dict()
        
        ξs = boundle.ξs
        βs = good_βs[stst]
        ider_map = data[stst]["ider_map"]
        
        # Measured flux related values
        for rath_ider in ["μ"; Rd.msd_mets]
            
            # ider map
            model_ider = rath_ider == "μ" ? 
                obj_ider : HG.exch_met_map[HG.mets_map[rath_ider]]
            
            ider_map[model_ider] = rath_ider
            ider_map[rath_ider] = model_ider
            
            for ξ in ξs
                # FBA
                data[stst][(ξ, :fba, rath_ider)] = (av = Ch.Utils.av(boundle, ξ, :fba, model_ider),)
                # bounds
                model = boundle[ξ, :net]
                data[stst][(ξ, :bounds, rath_ider)] = Ch.Utils.bounds(model, model_ider)
                
                for β in βs
                    
                    # EP
                    data[stst][(ξ, β, :ep, rath_ider)] = (
                        av = Ch.Utils.av(boundle, ξ, β, :ep, model_ider),
                        va = Ch.Utils.va(boundle, ξ, β, :ep, model_ider),
                        μ = Ch.Utils.μ(boundle, ξ, β, :ep, model_ider),
                        σ = Ch.Utils.σ(boundle, ξ, β, :ep, model_ider)
                    )

                end
            end
            
        end # rath_ider for
        
        # Stoi error
        for ξ in boundle.ξs
            metnet = boundle[ξ, :net]
            epouts = boundle[ξ, good_βs[stst], :ep]
            errs = map(epouts) do epout 
                Ch.Utils.norm_abs_stoi_err(metnet, epout)
            end
            data[stst][(ξ, "mean_ep_norm_stoi_err")] = mean.(errs)
            data[stst][(ξ, "max_ep_norm_stoi_err")] = maximum.(errs)
        end
        
        
    end # boundles for
    
    # Save
    file = joinpath(tIG.MODEL_PROCESSED_DATA_DIR, splitext(basename(dat_file))[1] * "___extracted_data.json")
    JSON.write(file, JSON.json(data))
    println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
end


