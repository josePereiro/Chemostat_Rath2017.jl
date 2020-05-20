# DIRS
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, MODEL_NAME)
mkpath(MODEL_PROCESSED_DATA_DIR)
const MODEL_RAW_MAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "$(MODEL_NAME)_url.mat")
Base.include_dependency(MODEL_RAW_MAT_FILE)

# FILES
# csv
const METS_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "mets_map.csv")
Base.include_dependency(METS_MAP_FILE)
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.csv")
Base.include_dependency(EXCH_MET_MAP_FILE)
const BASE_INTAKE_INFO_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_intake_info.csv")
Base.include_dependency(BASE_INTAKE_INFO_FILE)

# jls
const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_model.jls")
const FVA_PP_BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_preprocessed_base_model.jls")