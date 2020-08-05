# DIRS
const MODEL_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Human1_Publication_Data_Scripts")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, "Human1", MODEL_NAME)
const MODEL_CACHE_DATA_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")
const MODEL_FIGURES_DATA_DIR = joinpath(FIGURES_DATA_DIR, "Human1", MODEL_NAME)

function _create_dirs()
    mkpath(MODEL_PROCESSED_DATA_DIR)
    mkpath(MODEL_CACHE_DATA_DIR)
    mkpath(MODEL_FIGURES_DATA_DIR)
end


# FILES
# mat 
const MODEL_RAW_MAT_FILE = joinpath(MODEL_RAW_DATA_DIR, "tINIT_GEMs/data/HumanGEM.mat")

# csv
const METS_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "mets_map.bson")
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.bson")
const NIKLAS_BIOMASS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "niklas_biomass.bson")
const BASE_INTAKE_INFO_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_intake_info.bson")
const BASE_READABLE_MET_IDS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "readable_met_ids.bson")
const HAM_MEDIUM_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "ham_medium.bson")
const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_model.bson")


# jls
const FVA_PP_BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_preprocessed_base_model.bson")