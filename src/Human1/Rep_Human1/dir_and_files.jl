# DIRS
const MODEL_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Human1_Publication_Data_Scripts")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, "Human1", NAME)
const MODEL_CACHE_DATA_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")
const MODEL_FIGURES_DATA_DIR = joinpath(FIGURES_DATA_DIR, "Human1", NAME)
const ECMODELS_DATA_DIR = joinpath(MODEL_RAW_DATA_DIR, "ec_GEMs/models")

function _create_dirs()
    mkpath(MODEL_PROCESSED_DATA_DIR)
    mkpath(MODEL_CACHE_DATA_DIR)
    mkpath(MODEL_FIGURES_DATA_DIR)
end

# FILES
const COMP_FVA_HG_INPUT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "compareFVA_humanGEM_input.bson")
const COMP_FVA_HG_OUTPUT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "compareFVA_humanGEM_output.bson")
const COMP_FVA_HG_EXTRACTED_DATA_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "compareFVA_humanGEM_extracted_data.bson")