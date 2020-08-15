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
quickactivate(@__DIR__, "Chemostat_Rath2017")

using SparseArrays

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, ecGEMs, Rep_Human1
# -

# TODO: create a general template by combining all availables ecModels in raw data
input_dat = wload(Rep_Human1.COMP_FVA_HG_INPUT_FILE);
ec_model = Ch.Utils.uncompress_model(input_dat["ec_model"]);

ec_template = ec_model;

# Saving
file = ecGEMs.MODEL_EC_TEMPLATE_FILE
ec_template = Ch.Utils.compress_model(ec_template);
tagsave(file, Dict(DATA_KEY => ec_template))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
