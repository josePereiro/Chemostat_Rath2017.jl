const MODEL_NAME = "MODEL1105100000"
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, MODEL_NAME)
const MODEL_RAW_MAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "$(MODEL_NAME)_url.mat")

# Cheking for mat file
if !isfile(MODEL_RAW_MAT_FILE)
    error("$(MODEL_RAW_MAT_FILE) not found, you must run src/MODEL1105100000/1_make_mat_file.py fisrt!!!")
end