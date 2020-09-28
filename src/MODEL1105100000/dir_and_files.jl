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

# csv
const METS_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "mets_map.bson")
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.bson")
const BASE_INTAKE_INFO_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_intake_info.bson")
const NIKLAS_BIOMASS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "niklas_biomass.bson")

# jls
const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_model.bson")
const FVA_PP_BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "fva_preprocessed_base_model.bson")
const MAXENT_FBA_EB_BUNDLES_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "maxent_fba_ep_bundle.bson")