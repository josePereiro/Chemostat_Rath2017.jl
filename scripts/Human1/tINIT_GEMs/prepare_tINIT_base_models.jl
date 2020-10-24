# -*- coding: utf-8 -*-
# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import JSON
import MAT
using SparseArrays

import Chemostat
const Ch = Chemostat

import Chemostat_Rath2017: DATA_KEY, RathData, Human1
import Chemostat_Rath2017.Human1: HumanGEM, tINIT_GEMs
const HG = HumanGEM
const tIG = tINIT_GEMs
const Rd = RathData;
# -

# ---
# ## Description
# ---

# This script will load all the raw iINIt models and prepare then for modeling. See comments in the code for details

# ---
# ## Load brain models
# ---

tINIT_brain_models = wload(tIG.tINIT_RAW_BRAIN_MODELS_FILE)[DATA_KEY];

# TODO: fill with explanations
tINIT_base_models = Dict()
for (dat_id, tINIT_dat) in tINIT_brain_models
    println("\n ----------------- Test Processing $dat_id -----------------\n")
    base_model = tINIT_dat["metnet"]
    base_model = Ch.Utils.uncompress_model(base_model)
    mat_model = tINIT_dat["mat"]
    Ch.Utils.clamp_bounds!(base_model, Human1.MAX_BOUND, Human1.ZEROTH);
    
    # fix b dimention
    empty!(base_model.b) 
    push!(base_model.b, zeros(size(base_model, 1))...)

    base_model = Human1.delete_boundary_mets(base_model)

    exchs = Human1.prepare_extract_exchanges!(base_model)
    
    medium = HG.base_intake_info |> keys |> collect
    medium = filter((rxn) -> rxn in base_model.rxns, medium)
    Human1.open_rxns!(base_model, medium)
    
    Human1.apply_biomass!(base_model, Human1.BIOMASS_IDER, HG.niklas_biomass)
    base_model = Human1.set_atp_demand(base_model)
    
    Human1.try_fba(base_model, Human1.BIOMASS_IDER)
    
    fbaout = Human1.try_fba(base_model, Human1.BIOMASS_IDER);
    @assert fbaout.obj_val > Human1.ZEROTH
    
    base_model = Ch.Utils.compress_model(base_model)
    tINIT_base_models[dat_id] = Dict()
    tINIT_base_models[dat_id]["metnet"] = base_model
    tINIT_base_models[dat_id]["mat"] = mat_model
end

# Saving
file = tIG.tINIT_BASE_BRAIN_MODELS_FILE
tagsave(file, Dict(DATA_KEY => tINIT_base_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


