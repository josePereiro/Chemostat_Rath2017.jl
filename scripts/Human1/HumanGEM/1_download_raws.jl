## --------------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## --------------------------------------------------------------------
import Chemostat_Rath2017
const ChR = Chemostat_Rath2017
const HG = ChR.Human1.HumanGEM

## --------------------------------------------------------------------
# const zenodo_download_link = "https://zenodo.org/record/3583004/files/Human1_Publication_Data_Scripts.zip?download=1"
zenodo_download_link = "https://zenodo.org/record/3583004/files/Human1_Publication_Data_Scripts.zip"

cd(HG.MODEL_RAW_DATA_DIR |> dirname)
println("Now at: ", pwd())
zip_file = "Human1_Publication_Data_Scripts.zip"
if isdir(HG.MODEL_RAW_DATA_DIR)
    println(relpath(HG.MODEL_RAW_DATA_DIR), " already exist, to force a re-unzip delete the folder")
elseif isfile(zip_file)
    println(relpath(zip_file), " already exist, to force a re-download delete it")
    run(`unzip $zip_file`)
else
    println("Downloading using curl ...")
    run(`curl $(zenodo_download_link) --output $(zip_file)`)
    run(`unzip $zip_file`)
end
!isdir(HG.MODEL_RAW_DATA_DIR) && 
    error("$(HG.MODEL_RAW_DATA_DIR) not found after download && unzip!!!")
println(relpath(HG.MODEL_RAW_DATA_DIR), " ready!!!")
flush(stdout)
flush(stderr)