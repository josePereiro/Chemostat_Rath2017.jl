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
#=
From https://doi.org/10.5281/zenodo.3577466
ecGEMs/README.md
## Flux variability analysis
Flux variability analysis (FVA) (corresponding to the results presented in Fig. 5B) can 
be run using the `comparativeFVA_humanModels.m` function in the `ComplementaryScripts/Simulation` 
subdirectory. Specify the name of the model (cell line) for which FVA is to be run; for example:

This scripts reproduce this analysis
=#

# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import MAT
using StatsBase

import Chemostat
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, HumanGEM, RathData, Rep_Human1, 
                            print_action, load_cached, save_cache, set_cache_dir,
                            delete_temp_caches, temp_cache_file
const HG = HumanGEM
const Rd = RathData
const RepH1 = Rep_Human1;
# -

# ### Loading Models
# ---

CELL_LINE_NAME = "HOP62"
CELL_LINE_DIR = joinpath(RepH1.ECMODELS_DATA_DIR, CELL_LINE_NAME)
@assert isdir(CELL_LINE_DIR)
tINIT_MODEL_FILE = joinpath(CELL_LINE_DIR, CELL_LINE_NAME * ".mat");
@assert isfile(tINIT_MODEL_FILE)
EC_MODEL_FILE = joinpath(CELL_LINE_DIR, "ecModel_batch" * ".mat");
@assert isfile(EC_MODEL_FILE)

orig_dat = MAT.matread(tINIT_MODEL_FILE)["HOP62"];
orig_model = Ch.Utils.MetNet(orig_dat);
println("Orig model: ", size(orig_model))

ec_dat = MAT.matread(EC_MODEL_FILE)["ecModel_batch"]
ec_model = Ch.Utils.MetNet(ec_dat);
println("EC model: ", size(ec_model))

# ### Preparing Models
# ---

# globals
obj_ider = "biomass_human"
prot_pool_exchange = "prot_pool_exchange"
zeroth = 1e-8; # a minimum threshold to not be consider zero

# Remove boundary metabolites
println("\nRemoving boundary metabolites")
for var in [:orig_model, :ec_model]
    println("\tModel: ", var)
    model = eval(var)
    println("\tBefore size: ", size(model))
    to_remove = findall(model.mets) do met
        endswith(met, "x") # boundary metabolites ids ends in 'x'
    end
    println("\tTo delete: ", length(to_remove))
    model = Ch.Utils.del_met(model, to_remove)
    println("\tAfter size: ", size(model))
    eval(:(global $var = $model))
    println()
end

# +
mediaComps =["glucose", "arginine", "histidine", "lysine", "methionine", "phenylalanine", "tryptophan", 
            "tyrosine", "alanine", "glycine", "serine", "threonine", "aspartate", "glutamate", "asparagine", 
            "glutamine", "isoleucine", "leucine", "proline", "valine", "cysteine", "thiamin", "hypoxanthine", 
            "folate", "biotin", "pantothenate", "choline", "inositol", "nicotinamide", "pyridoxine", 
            "riboflavin", "thymidine", "aquacob(III)alamin", "lipoic acid", "sulfate", "linoleate", 
            "linolenate", "O2", "H2O", "retinoate", "Fe2+", "Pi", "alpha-tocopherol", "gamma-tocopherol"];
mediaComps .*= "[s]"; # External met readable names ends in '[s]'

# Get echange reactions (compatible with original HumanGEM ids)
mediumExcIds = [HG.exch_met_map[HG.readable_met_ids_map[comp]] for comp in mediaComps];

# +
# Open all medium exchanges to (-1000, 1000) and close all the other exchanges intake lb = 0
println("\nProcessing exchanges")

# Orig model
println("\tOrig model")
# Get all the exchanges 
# Exchange reactions are defined as reactions which involve only products or only reactants
exchs = Ch.Utils.exchanges(orig_model)
println("\tExchs: ", length(exchs))
foreach(exchs) do rxn
    Ch.Utils.lb!(orig_model, rxn, 0.0);
end
foreach(mediumExcIds) do rxn
    if rxn in orig_model.rxns
        Ch.Utils.lb!(orig_model, rxn, -1000);
        Ch.Utils.ub!(orig_model, rxn, 1000);
    end
end

# report unused metabolites
foreach(mediumExcIds) do rxn
    bounds = Ch.Utils.bounds(orig_model, rxn)
    if bounds != (-1000, 1000)
        @warn(rxn, " not used as medium component!!!")
    end
end

println()

# +
# Ec model
# The Ec model is in irreversible format, meaning that a reaction A <-> B 
# is splitted in two reactions A -> B, B -> A
println("\tEc model")
exchs = Ch.Utils.exchanges(ec_model)
println("\tExchs: ", length(exchs))
for rxni in exchs
    rxn = ec_model.rxns[rxni]
    
    # Exclude protein pool exchange
    rxn == prot_pool_exchange && continue
    
    # Differentiate between uptakes and production reactions
    if endswith(rxn, "_REV") # is uptake
        # close all uptakes
        Ch.Utils.lb!(ec_model, rxn, 0.0);
        Ch.Utils.ub!(ec_model, rxn, 0.0);
        
        # if it is in the medium open it
        real_id = replace(rxn, "_REV" => "")
        if real_id in mediumExcIds
            Ch.Utils.lb!(ec_model, rxn, 0.0);
            Ch.Utils.ub!(ec_model, rxn, 1000.0);
        end
            
    else # is production
        # Open all production reactions
        Ch.Utils.lb!(ec_model, rxn, 0.0);
        Ch.Utils.ub!(ec_model, rxn, 1000);
    end
end

# report unused metabolites
foreach(mediumExcIds) do rxn
    model_rxn = rxn * "_REV"
    bounds = Ch.Utils.bounds(ec_model, model_rxn)
    if bounds != (0.0, 1000)
        @warn(rxn, " not used as medium component!!!")
    end
end

println()
# -

# Gets the optimal value for ecirrevModel and fixes the objective value to
# this for both models
println("\nFixing objective ($obj_ider)")
obj_val = Ch.LP.fba(ec_model, obj_ider).obj_val
println("\tobj_val: ", obj_val)
for var in [:orig_model, :ec_model]
    model = eval(var)
    Ch.Utils.bounds!(model, obj_ider, obj_val - zeroth, obj_val + zeroth)
end

# TODO: implement parsimonious fba
# Get a parsimonious (in this implementation we do not make it parsimonious) 
# flux distribution for the ecModel (minimization of
# total protein usage)
obj_val = Ch.LP.fba(ec_model, prot_pool_exchange; sense = 1.0).obj_val
Ch.Utils.bounds!(ec_model, prot_pool_exchange, obj_val - zeroth, obj_val + zeroth)

println("\n Check if models are feasible")
for var in [:orig_model, :ec_model]
    println(var)
    model = eval(var)
    obj_val = Ch.LP.fba(model, obj_ider).obj_val;
    if obj_val < zeroth
        println("\tWARNING: Constrained $var is unfeasible")
    else
        println("\tConstrained $var is feasible, objval: ", obj_val)
    end
end

## Saving input
println("\nSaving")
file = RepH1.COMP_FVA_HG_INPUT_FILE
tagsave(file, Dict( 
    "orig_model" => Ch.Utils.compress_model(orig_model),
    "ec_model" => Ch.Utils.compress_model(ec_model)
    )
)
println(relpath(file), " created, size: ", filesize(file), " bytes.")

