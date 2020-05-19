# Here I group all the computation related with the MODEL1105100000 model from 
# Shlomi et al., Genome-Scale Metabolic Modeling Elucidates the Role of 
# Proliferative Adaptation in Causing the Warburg Effect.")
# TODO search model download link

module MODEL1105100000
    using ..Chemostat_Rath2017: PROJ_ROOT, PROCESSED_DATA_DIR
    import CSV
    import DataFrames: DataFrame

    include("1_meta.jl")
    include("2_mets_map.jl")
    # include("2_rxns_map.py")
end

