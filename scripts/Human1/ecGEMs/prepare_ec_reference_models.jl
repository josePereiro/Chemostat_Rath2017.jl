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

import MAT
using SparseArrays
using Test
using Distributions

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Human1
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
# -

# ---
# ### Collect files

# Pick all ec Models
data_dir = joinpath(ecG.MODEL_RAW_DATA_DIR, "ec_GEMs/models") # TODO: package this
ec_model_files = []
for (root, _, files) in walkdir(data_dir)
    for file in files
        file = joinpath(root, file)
        if basename(file) == "ecModel_batch.mat"
            push!(ec_model_files, file)
        end
    end
end

# ---
# ### Processing

ec_reference_models = []
for (i, file) in ec_model_files |> enumerate
    model = Ch.Utils.read_mat(file);
    Ch.Utils.clamp_bounds!(model, Human1.MAX_BOUND, Human1.ZEROTH);
    name = file |> dirname |> basename
    println("\n\n--------------- model $name [$i/$(length(ec_model_files))]---------------\n")
    
    # We delete the boundary metabolites, they are not required
    # bkwd and fwd splitted reactions are troublemakers for EP, but they 
    # are necessary to model enzymatic costs. So, we leave as least as possible 
    # only leaving in rev format the exchanges

    # We unified the exchanges (make them a unique rxn), and let them in a 
    # semi-open state (intake bloked, outtake open)   
    model = Human1.delete_boundary_mets(model);

    # We close all reactions marked as "Exchanges/boundary" and returns the exchanges
    exchs, bkwd_exchs, fwd_exchs = Human1.prepare_extract_exchanges!(model, Human1.MAX_BOUND);

    # We delate any backward defined exchange reaction, 
    # all reactions will be possible reversible
    model, exchs, bkwd_exchs, fwd_exchs = Human1.del_bkwd_exchs(model, bkwd_exchs);

    # Apply Hams medium
    for (rxn, c) in HG.base_intake_info
        !(rxn in model.rxns) && continue
        Ch.Utils.bounds!(model, rxn, -Human1.MAX_BOUND, Human1.MAX_BOUND);
    end

    # Open prot pool exchange
    Ch.Utils.bounds!(model, Human1.PROT_POOL_EXCHANGE, 0.0, Human1.MAX_BOUND);

    # Check that the model is feasible
    fbaout = Human1.try_fba(model, Human1.OBJ_IDER);
    @assert fbaout.obj_val > Human1.ZEROTH
    
    model = Ch.Utils.compress_model(model);
    push!(ec_reference_models, model)
end

# saving
file = ecG.EC_BASE_REFERENCE_MODELS
tagsave(file, Dict(DATA_KEY => ec_reference_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
