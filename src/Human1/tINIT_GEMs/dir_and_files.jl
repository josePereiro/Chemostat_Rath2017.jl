# DIRS
const MODEL_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Human1_Publication_Data_Scripts")
const RAW_tINIT_OUTPUTS_DIR = joinpath(MODEL_RAW_DATA_DIR, "tINIT_GEMs", "run_tINIT_outputs")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, MODEL_NAME)
mkpath(MODEL_PROCESSED_DATA_DIR)
const MODEL_CACHE_DATA_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")
mkpath(MODEL_CACHE_DATA_DIR)
const MODEL_FIGURES_DATA_DIR = joinpath(FIGURES_DATA_DIR, MODEL_NAME)
mkpath(MODEL_FIGURES_DATA_DIR)

# FILES
# mat 
const MODEL_RAW_MAT_FILE = joinpath(MODEL_RAW_DATA_DIR, "HumanGEM.mat")
Base.include_dependency(MODEL_RAW_MAT_FILE)

# # jls
# const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_model.jls")
# const FVA_PP_BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_preprocessed_base_model.jls")