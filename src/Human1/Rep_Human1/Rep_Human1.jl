"""
    This module assist the scripts that reproduce the work on
    Robinson, Jonathan L., Pınar Kocabaş, Hao Wang, Pierre-Etienne Cholley, Daniel Cook, Avlant Nilsson, Mihail Anton, et al. “An Atlas of Human Metabolism.” 
    Science Signaling 13, no. 624 (March 24, 2020). https://doi.org/10.1126/scisignal.aaz1482.
"""
module Rep_Human1
    import ...Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR, FIGURES_DATA_DIR, RathData

    include("meta.jl")
    include("dir_and_files.jl")

    function __init__()
        _create_dirs()
    end
end