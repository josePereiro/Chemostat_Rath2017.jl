# DIRS
const MODEL_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, MODEL_NAME)
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

# csv
const METS_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "mets_map.csv")
Base.include_dependency(METS_MAP_FILE)
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.csv")
Base.include_dependency(EXCH_MET_MAP_FILE)
