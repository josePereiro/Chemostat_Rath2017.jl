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
@quickactivate "Chemostat_Rath2017"

import MAT
using SparseArrays
using Test
using Distributions
using ProgressMeter

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017
import Chemostat_Rath2017: DATA_KEY, Human1
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
# -

# ---
# ## Description
# This script take the given ec_model as a template to produce a new one from the orig_model

# load orig input models
# tINIT_brain_models = wload(tIG.tINIT_RAW_BRAIN_MODELS_FILE)[DATA_KEY];
ec_refdata_set = wload(ecG.EC_REFERENCE_DATA)[DATA_KEY] |> values |> collect;

# +
# rxn_ecmap = ec_refdata["HOP92"]["ecmap"];
# -

brain_models_dat = wload(tIG.tINIT_BASE_BRAIN_MODELS_FILE)[DATA_KEY];

model = brain_models_dat["GTEx-brain"]["metnet"]
model = Ch.Utils.uncompress_model(model);

ec_model = Human1.build_ecModel(model, ec_refdata_set);
@time Human1.try_fba(ec_model, Human1.OBJ_IDER);













 























test_model, allowed_protless, kin_stoi, draw_stoi = build_ecModel(model, [ec_refdata["HOP92"]]);

test_model = deepcopy(ec_model);
# Human1.try_fba(test_model, Human1.OBJ_IDER);
test_model = merge_protless!(test_model, allowed_protless, kin_stoi, draw_stoi);
@time Human1.try_fba(test_model, Human1.OBJ_IDER);



# +
# const expanded_model = Human1.expanded_model
# const ECMAP_KEY = Human1.ECMAP_KEY
# const PROTLESS_KEY = Human1.PROTLESS_KEY
# const PROT_STOIS_KEY = Human1.PROT_STOIS_KEY
# const EXTRAS_KEY = Human1.EXTRAS_KEY
# const merge_rxndata! = Human1.merge_rxndata!
# const RxnData = Human1.RxnData
# const clear_rxndata! = Human1.clear_rxndata!
# const compated_model = Human1.compated_model
# const collect_protless = Human1.collect_protless
# const EMPTY_SPOT = Human1.EMPTY_SPOT
# const REV_SUFFIX = Human1.REV_SUFFIX
# # const free_spot = Human1.free_spot
# -

Ch.Utils.summary(test_model, Human1.PROT_POOL_EXCHANGE)
Ch.Utils.summary(ec_model, Human1.PROT_POOL_EXCHANGE)

Ch.Utils.summary(test_model, "draw_prot_A1L3X0")
Ch.Utils.summary(ec_model, "draw_prot_A1L3X0")
Ch.Utils.summary(test_model, "draw_prot_11DOCRTSTRNte")

Ch.Utils.summary(ec_model, "EX_lvstacid[e]")

exchs = Ch.Utils.exchanges(ec_model);

for rxn in exchs
    ec_model.subSystems[rxn] |> println
end



all_protless = ec_model.rxns[collect_protless(ec_model)]
disallowed_protless = setdiff(all_protless, allowed_protless) |> unique



allowed_protless

Human1.try_fba(ec_model, Human1.OBJ_IDER);

# ---
# ## Build ec template

ec_brain_models = Dict()
tINIT_brain_models = [tINIT_brain_models["GTEx-brain"]] # Test
for (model_id, dat) in tINIT_brain_models
    
    println("\n\n ------------- Processing $model_id ------------- \n\n")

    model = dat["metnet"]
    println("Orig model size: ", size(model))
    
    gd = Human1.EcMGDC(model, ec_template);
    Human1.collect_rxns_data!(gd);
    Human1.collect_mets_data!(gd);
    Human1.collect_draw_rxns!(gd);
    Human1.add_prot_pool_exchange!(gd);
    Human1.make_all_unique!(gd);
    flush(stdout)
            
    ec_model = Human1.build_new_model(gd);
    ec_model = Ch.Utils.compress_model(ec_model)
    println("Ec model size: ", size(ec_model))
    ec_brain_models[model_id] = ec_model
    
    flush(stdout)
end

# +
# TODO: Ramake the ec inclusion workflow to process one reaction at a time, and possible 
# from several ec reference models, this will make possible to check fba at any time.

# +
# # saving
# file = ecG.EC_BRAIN_RAW_MODELS_FILE
# tagsave(file, Dict(DATA_KEY => ec_brain_models))
# println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
# -


