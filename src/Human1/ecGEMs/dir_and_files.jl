# DIRS
const MODEL_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Human1_Publication_Data_Scripts")
const RAW_tINIT_OUTPUTS_DIR = joinpath(MODEL_RAW_DATA_DIR, "tINIT_GEMs", "run_tINIT_outputs")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, "Human1" , PROJ_NAME)
const MODEL_CACHE_DATA_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")
const MODEL_FIGURES_DATA_DIR = joinpath(FIGURES_DATA_DIR, PROJ_NAME)
const EC_RAW_MODELS_DIR = joinpath(MODEL_RAW_DATA_DIR, "ec_GEMs/models")
function _create_dirs()
    for dir in [MODEL_PROCESSED_DATA_DIR, MODEL_CACHE_DATA_DIR, MODEL_FIGURES_DATA_DIR]
        if !isdir(dir)
            mkpath(dir)
        end
    end
end

# FILES
# mat 
const EC_BRAIN_BASE_MODELS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "Human1_ec_base_brain_models.bson")
const EC_REFERENCE_DATA = joinpath(MODEL_PROCESSED_DATA_DIR, "ec_reference_data.bson")
const FVA_PP_BASE_MODELS = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_pp_models.bson")
const MAXENT_FBA_EB_BOUNDLES_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "maxent_fba_ep_bundles.bson")
const EXTRACTED_DATA_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "extracted_data.bson")

