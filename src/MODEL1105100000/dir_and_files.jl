# DIRS
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, PROJ_NAME)
mkpath(MODEL_PROCESSED_DATA_DIR)
const MODEL_CACHE_DATA_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")
mkpath(MODEL_CACHE_DATA_DIR)
const MODEL_FIGURES_DATA_DIR = joinpath(FIGURES_DATA_DIR, PROJ_NAME)
mkpath(MODEL_FIGURES_DATA_DIR)


# FILES
# mat 
const MODEL_RAW_MAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "$(PROJ_NAME)_url.mat")
Base.include_dependency(MODEL_RAW_MAT_FILE)

# csv
const METS_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "mets_map.csv")
Base.include_dependency(METS_MAP_FILE)
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.csv")
Base.include_dependency(EXCH_MET_MAP_FILE)
const BASE_INTAKE_INFO_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_intake_info.csv")
Base.include_dependency(BASE_INTAKE_INFO_FILE)
const NIKLAS_BIOMASS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "niklas_biomass.csv")
Base.include_dependency(NIKLAS_BIOMASS_FILE)

# jls
const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_model.jls")
const FVA_PP_BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_preprocessed_base_model.jls")