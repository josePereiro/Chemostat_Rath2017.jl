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

using SparseArrays
import JSON

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017
import Chemostat_Rath2017: DATA_KEY, RathData, Human1
import Chemostat_Rath2017.Human1: Rep_Human1, ecGEMs, tINIT_GEMs, HumanGEM
const RepH1 = Rep_Human1;
const ecG = ecGEMs
const tIG = tINIT_GEMs;
const HG = HumanGEM
const Rd = RathData
# -

# ---
# ## Loading data

println("\nLoading brain models")
@time brain_models = wload(tIG.tINIT_BASE_BRAIN_MODELS_FILE)[DATA_KEY];

println("\nLoading ec reference data set")
@time ec_refdata_dict = wload(ecG.EC_REFERENCE_DATA)[DATA_KEY];
ec_refdata_set = ec_refdata_dict |> values |> collect;

# ---
# ## Generating EC Models

ec_models = Dict()
for (model_id, dat) in brain_models
    println("\n\n--------------- Generating $model_id ecModel ---------------")
    
    src_model = dat["metnet"]
    src_model = Ch.Utils.uncompress_model(src_model)
    println("src model")
    println("size: ", size(src_model))
    fbaout = Human1.try_fba(src_model, Human1.OBJ_IDER)
    @assert fbaout.obj_val > Human1.ZEROTH
    
    ec_model = Human1.build_ecModel(src_model, ec_refdata_set; add_protless = true);
    println("\nec model")
    println("size: ", size(ec_model))
    fbaout = Human1.try_fba(ec_model, Human1.OBJ_IDER)
    @assert fbaout.obj_val > Human1.ZEROTH
    Human1.print_ec_stats(ec_model)
    
    ec_models[model_id] = ec_model 
end

# ---
# ## Prepare base models

base_models = Dict()
for (model_id, ec_model) in ec_models
    
    println("\n\n--------------- Generating $model_id ec base model ---------------")
    
    base_model = Human1.set_atp_demand(ec_model)
#     base_model = deepcopy(model)
    Ch.Utils.clamp_bounds!(base_model, Human1.MAX_BOUND, Human1.ZEROTH);
    
    exchs = Human1.prepare_extract_exchanges!(base_model)

    medium = HG.base_intake_info |> keys |> collect
    medium = filter((rxn) -> rxn in base_model.rxns, medium)
    Human1.open_rxns!(base_model, medium)

    Human1.apply_biomass!(base_model, Human1.OBJ_IDER, HG.niklas_biomass)

    fbaout = Human1.try_fba(base_model, Human1.OBJ_IDER);
    @assert fbaout.obj_val > Human1.ZEROTH
    
    base_models[model_id] = Ch.Utils.compress_model(base_model)
end

# saving
# TODO: save all as compressed dicts
# file = ecG.EC_BRAIN_BASE_MODELS_FILE
# tagsave(file, Dict(DATA_KEY => base_models))
# println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


