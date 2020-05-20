# Here I group all the computation related with the MODEL1105100000 model from 
# Shlomi et al., Genome-Scale Metabolic Modeling Elucidates the Role of 
# Proliferative Adaptation in Causing the Warburg Effect.")
# TODO search model download link

module MODEL1105100000
    import ..Chemostat_Rath2017: PROJ_ROOT, PROCESSED_DATA_DIR
    # import .RathData
    # Rd = RathData
    import CSV
    # import MAT
    import Serialization: deserialize
    import DataFrames: DataFrame
    # import Chemostat
    # Ch = Chemostat
    

    include("meta.jl")
    include("dir_and_files.jl")
    include("load_data.jl")
end