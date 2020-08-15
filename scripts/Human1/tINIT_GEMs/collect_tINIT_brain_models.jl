# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl
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

import MAT

import Chemostat.Utils: MetNet, compress_model
import Chemostat_Rath2017: DATA_KEY, tINIT_GEMs
const tIG = tINIT_GEMs;
# -

# ---
# ### HumanGEM tINIT_output files

# Those are all the tINIT_output data related with HumanGEM
mat_files = Dict()
mat_files["GTEx"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "GTEx/tINIT_GTEx_outputs.mat");
# mat_files["Hart2015"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "Hart2015/tINIT_Hart2015_HumanGEM_outputs.mat");
mat_files["TCGA"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "TCGA/tINIT_TCGA_outputs.mat")
# mat_files["DepMap_1_1"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_1_1.mat")
# mat_files["DepMap_1_2"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_1_2.mat")
# mat_files["DepMap_2_1"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_2_1.mat")
# mat_files["DepMap_2_2"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_2_2.mat")
# mat_files["DepMap_3_1"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_3_1.mat")
# mat_files["DepMap_3_2"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_3_2.mat")
# mat_files["DepMap_4_1"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_4_1.mat")
# mat_files["DepMap_4_2"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_4_2.mat")
# mat_files["DepMap_5_1"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_5_1.mat")
# mat_files["DepMap_5_2"] = joinpath(tIG.RAW_tINIT_OUTPUTS_DIR, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_5_2.mat");
dataset_ids = mat_files |> keys |> collect |> sort;

# ---
# ### Collecting models

# ids and names taken from data/raw/Human1_Publication_Data_Scripts/tINIT_GEMs/figures/model_id_type_mapping.txt
healty_models_id = ["brain", "GBM NT"];
cancer_models_id = ["brain-cancer", "GBM TP", "GBM TR", "LGG TP", "LGG TR"];

brain_models = Dict()
for dataset_id in dataset_ids
    
    println("Data set id: ", dataset_id)
    
    file = mat_files[dataset_id]
    println("\tFile: ", relpath(file))
    
    dat = MAT.matread(file)["INIT_output"]
    mat_models = dat["model"]
    ids = get(dat, "id", get(dat, "tissues", nothing))
    isnothing(ids) && error("ids not found!!!")
    
    for (model_id, mat_model) in zip(ids, mat_models)
        if model_id in healty_models_id || model_id in cancer_models_id
            model_key = "$dataset_id-$model_id"
            println("\tModel found: ", model_key); flush(stdout)
            
            brain_models[model_key] =  Dict()

            metnet = MetNet(mat_model)
            metnet = compress_model(metnet)
            brain_models[model_key]["metnet"] = metnet

            # lite up mat_model
            foreach(["S", "b", "lb", "ub", "c"]) do k
                delete!(mat_model, k)
            end
            brain_models[model_key]["mat"] = mat_model
            
        end
    end
    println()
end

# saving
file = tIG.tINIT_BRAIN_MODELS_FILE
tagsave(file, Dict(DATA_KEY => brain_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
