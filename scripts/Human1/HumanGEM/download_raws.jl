# +
import CSV
import DataFrames: DataFrame
import MAT

import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
# -

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

zenodo_download_link = "https://zenodo.org/record/3583004/files/Human1_Publication_Data_Scripts.zip?download=1"

if isdir(HG.MODEL_RAW_DATA_DIR)
    println(relpath(HG.MODEL_RAW_DATA_DIR), " already exist, to force a re-download delete the folder")
else
    zip_file = joinpath(HG.MODEL_RAW_DATA_DIR |> dirname, "Human1_Publication_Data_Scripts.zip")
    download(zenodo_download_link, zip_file)
    run(`unzip $zip_file $(HG.MODEL_RAW_DATA_DIR)`)
    !isdir(HG.MODEL_RAW_DATA_DIR) && error("$(HG.MODEL_RAW_DATA_DIR) not found after download!!!")
    println(relpath(HG.MODEL_RAW_DATA_DIR), " downloaded")
end
flush(stdout)
flush(stderr)