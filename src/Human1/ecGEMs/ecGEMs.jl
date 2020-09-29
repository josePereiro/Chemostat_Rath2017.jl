# Raw data
# The Human1 GEM-PRO, as well as all cell-specific GEMs and
# enzyme-constrained ecGEMs described in this study, and the datasets and scripts necessary
# to reproduce the models, analyses, and figures are available on Zenodo (https://doi.org/10.5281/zenodo.3577466). 

module ecGEMs
    import ...Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR, FIGURES_DATA_DIR, RathData
    
    include("const.jl")
    include("dir_and_files.jl")

    function __init__()
        _create_dirs()
    end
end