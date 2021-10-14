## --------------------------------------------------------------------
using ProjAssistant
quickactivate(@__DIR__, "Chemostat_Rath2017")

## --------------------------------------------------------------------
import Chemostat_Rath2017
const ChR = Chemostat_Rath2017
const HG = ChR.Human1

## --------------------------------------------------------------------
const ZENODO_DOWNLOAD_LINK = "https://zenodo.org/record/3583004/files/Human1_Publication_Data_Scripts.zip"

cd(rawdir(HG))
println("Now at: ", pwd())

const DEST_DIT = rawdir(HG, "Human1_Publication_Data_Scripts")
const ZIP_FILE_NAME = "Human1_Publication_Data_Scripts.zip"
if isdir(DEST_DIT)
    println(relpath(DEST_DIT), " already exist, to force a re-unzip delete the folder")
elseif isfile(ZIP_FILE_NAME)
    println(relpath(ZIP_FILE_NAME), " already exist, to force a re-download delete it")
    run(`unzip $ZIP_FILE_NAME`)
else
    println("Downloading using curl ...")
    run(`curl $(ZENODO_DOWNLOAD_LINK) --output $(ZIP_FILE_NAME)`)
    run(`unzip $(ZIP_FILE_NAME)`)
end
!isdir(DEST_DIT) && 
    error("$(DEST_DIT) not found after download && unzip!!!")
println(relpath(DEST_DIT), " ready!!!")
flush(stdout)
flush(stderr)