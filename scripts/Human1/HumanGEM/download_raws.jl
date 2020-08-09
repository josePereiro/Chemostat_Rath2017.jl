# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import CSV
import DataFrames: DataFrame
import MAT

import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
# -

zenodo_download_link = "https://zenodo.org/record/3583004/files/Human1_Publication_Data_Scripts.zip?download=1"

cd(HG.MODEL_RAW_DATA_DIR |> dirname)
zip_file = "Human1_Publication_Data_Scripts.zip"
if isdir(HG.MODEL_RAW_DATA_DIR)
    println(relpath(HG.MODEL_RAW_DATA_DIR), " already exist, to force a re-unzip delete the folder")
elseif isfile(zip_file)
    println(relpath(zip_file), " already exist, to force a re-download delete it")
    run(`unzip $zip_file`)
else
    download(zenodo_download_link, zip_file)
    run(`unzip $zip_file`)
end
!isdir(HG.MODEL_RAW_DATA_DIR) && error("$(HG.MODEL_RAW_DATA_DIR) not found after download-unzip!!!")
println(relpath(HG.MODEL_RAW_DATA_DIR), " ready!!!")
flush(stdout)
flush(stderr)