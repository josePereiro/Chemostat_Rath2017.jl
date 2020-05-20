# DIRS
RATH_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "rath2017___data")
mkpath(RATH_RAW_DATA_DIR)    
RATH_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(RATH_RAW_DATA_DIR))
mkpath(RATH_PROCESSED_DATA_DIR)    

# FILES
RATH_STDM_ORIG_FILE = joinpath(RATH_RAW_DATA_DIR, "rath2017___42_MAX_UB_standard_medium.tsv")
RATH_STDM_CONV_FILE = joinpath(RATH_PROCESSED_DATA_DIR, basename(RATH_STDM_ORIG_FILE))
Base.include_dependency(RATH_STDM_CONV_FILE)

RATH_CONT_CUL_DATA_FILE_SUFFIX = "rath2017___cont_exp";
RATH_CONT_CUL_DATA_ORIG_FILES = Dict()
RATH_CONT_CUL_DATA_CONV_FILES = Dict()
for exp in exps
    filename = "$(RATH_CONT_CUL_DATA_FILE_SUFFIX)_$(exp).tsv"
    RATH_CONT_CUL_DATA_ORIG_FILES[exp] = joinpath(RATH_RAW_DATA_DIR, filename)
    RATH_CONT_CUL_DATA_CONV_FILES[exp] = joinpath(RATH_PROCESSED_DATA_DIR, filename)
    Base.include_dependency(RATH_CONT_CUL_DATA_CONV_FILES[exp])
end

RATH_MAX_FLUX_ORIG_FILE = joinpath(RATH_RAW_DATA_DIR, "rath2017___max_invitro_fluxs.tsv")
RATH_MAX_FLUX_CONV_FILE = joinpath(RATH_PROCESSED_DATA_DIR, basename(RATH_MAX_FLUX_ORIG_FILE))
Base.include_dependency(RATH_MAX_FLUX_CONV_FILE)