# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
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
using DataFrames
using CSV

import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
H1 = Chemostat_Rath2017.Human1

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# +
## Ham's F-12 
## mediam concentration from from https://www.labome.com/method/Cell-Culture-Media-A-Review.html
## but mediam componets from Human1 Zenodo https://zenodo.org/record/3583004#.Xu8KdmpKiCh
## task (Growth in Ham's medium)
## repo TODO add link and explanation
## All concentrations in mM
Inf_ = 9999 # taked as infinite concentration
Ham_medium = Dict()

# Nutrients
Ham_medium["glucose[s]"] = 1e1
Ham_medium["lipoic acid[s]"] = 1.0e-3
Ham_medium["linoleate[s]"] = 3.0e-4
Ham_medium["hypoxanthine[s]"] = 3.0e-2
Ham_medium["thymidine[s]"] = 3.0e-3

# Essentials
Ham_medium["isoleucine[s]"] = 3.0e-2
Ham_medium["leucine[s]"] = 9.9e-2
Ham_medium["lysine[s]"] = 2.0e-1
Ham_medium["methionine[s]"] = 3.0e-2
Ham_medium["phenylalanine[s]"] = 3.0e-2
Ham_medium["threonine[s]"] = 1.0e-1
Ham_medium["tryptophan[s]"] = 9.8e-3
Ham_medium["valine[s]"] = 1.0e-1

# Semi-esentials
Ham_medium["cysteine[s]"] = 2.3e-1
Ham_medium["tyrosine[s]"] = 3.0e-2
Ham_medium["arginine[s]"] = 1e0
Ham_medium["histidine[s]"] = 1.0e-1

# Non-essentials
Ham_medium["alanine[s]"] = 1.0e-1
Ham_medium["asparagine[s]"] = 8.8e-2
Ham_medium["aspartate[s]"] = 1.0e-1
Ham_medium["glutamine[s]"] = 1.0e0
Ham_medium["glutamate[s]"] = 1.0e-1 
Ham_medium["glycine[s]"] = 1.0e-1
Ham_medium["proline[s]"] = 3.0e-1
Ham_medium["serine[s]"] = 1.0e-1


# Vitaminoids
Ham_medium["inositol[s]"] = 1.0e-10
Ham_medium["choline[s]"] = 1.0e-1

# Vitamins
Ham_medium["folate[s]"] = 2.0e-3
Ham_medium["pyridoxine[s]"] = 3.0e-3
Ham_medium["nicotinamide[s]"] = 3.0e-4 # Niacinamid in source
Ham_medium["pantothenate[s]"] = 1.0e-1
Ham_medium["biotin[s]"] = 3.0e-5
Ham_medium["riboflavin[s]"] = 1.0e-4
Ham_medium["thiamin[s]"] = 1.0e-3


# Supplements
Ham_medium["sulfate[s]"] = Inf_ 
Ham_medium["Fe2+[s]"] = Inf_ 
Ham_medium["Pi[s]"] = Inf_
Ham_medium["O2[s]"] = Inf_
Ham_medium["H2O[s]"] = Inf_


# Not found in concentration source
# TODO search real conc for this mets
Ham_medium["alpha-tocopherol[s]"] = 1e-2
Ham_medium["gamma-tocopherol[s]"] = 1e-2
Ham_medium["aquacob(III)alamin[s]"] = 1e-2
Ham_medium["retinoate[s]"] = 1e-3
Ham_medium["linolenate[s]"] = 0.0e0 # Not essential
Ham_medium["cysteine[s]"] = 0.0e0 # Not essential


# Saving
df = DataFrame(collect.([keys(Ham_medium), values(Ham_medium)]));
CSV.write(H1.HAM_MEDIUM_FILE, df)
println("created $(relpath(H1.HAM_MEDIUM_FILE))")
# -

# This names where taken from the task 'Growth on Ham's media' in the tINIT_GEMs/metabolic_tasks folder
Ham_medium_ids = ["arginine[s]", "histidine[s]", "lysine[s]", "methionine[s]", "phenylalanine[s]", 
    "tryptophan[s]", "tyrosine[s]", "alanine[s]", "glycine[s]", "serine[s]", "threonine[s]", "aspartate[s]", 
    "glutamate[s]", "asparagine[s]", "glutamine[s]", "isoleucine[s]", "leucine[s]", "proline[s]", "valine[s]", 
    "cysteine[s]", "thiamin[s]", "hypoxanthine[s]", "folate[s]", "biotin[s]", "pantothenate[s]", "choline[s]", 
    "inositol[s]", "nicotinamide[s]", "pyridoxine[s]", "riboflavin[s]", "thymidine[s]", "aquacob(III)alamin[s]",
    "lipoic acid[s]", "glucose[s]", "sulfate[s]", "linoleate[s]", "linolenate[s]", "O2[s]", "H2O[s]", 
    "retinoate[s]", "Fe2+[s]", "Pi[s]", "alpha-tocopherol[s]", "gamma-tocopherol[s]"];


