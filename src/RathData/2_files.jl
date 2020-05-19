RATH_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "rath2017___data")
if !isdir(RATH_RAW_DATA_DIR)
    error("$(RATH_RAW_DATA_DIR) not found!!!")
end
RATH_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(RATH_RAW_DATA_DIR))
if !isdir(RATH_PROCESSED_DATA_DIR)
    mkpath(RATH_PROCESSED_DATA_DIR)
    println("created $(RATH_PROCESSED_DATA_DIR)!!!")
end

RATH_STDM_ORIG_FILE = joinpath(RATH_RAW_DATA_DIR, "rath2017___42_MAX_UB_standard_medium.tsv")
if !isfile(RATH_STDM_ORIG_FILE)
    error("$(RATH_STDM_ORIG_FILE) not found!!!")
end
RATH_STDM_CONV_FILE = joinpath(RATH_PROCESSED_DATA_DIR, basename(RATH_STDM_ORIG_FILE))

RATH_CONT_CUL_DATA_FILE_SUFFIX = "rath2017___cont_exp";