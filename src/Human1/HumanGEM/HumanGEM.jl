# Raw data
# The Human1 GEM-PRO, as well as all cell-specific GEMs and
# enzyme-constrained ecGEMs described in this study, and the datasets and scripts necessary
# to reproduce the models, analyses, and figures are available on Zenodo (https://doi.org/10.5281/zenodo.3577466). 


module HumanGEM
    import BSON
    import DrWatson: wload
    import ...Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR, FIGURES_DATA_DIR, RathData, DATA_KEY
    import CSV
    import Serialization: deserialize
    import DataFrames: DataFrame

    include("meta.jl")
    include("dir_and_files.jl")
    include("load_data.jl")
    include("base_intake_info.jl")

    function __init__()
        _load_all_data()
        _create_dirs()
    end
end