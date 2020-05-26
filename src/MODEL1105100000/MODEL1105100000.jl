# Here I group all the computation related with the MODEL1105100000 model from 
# Shlomi et al., Genome-Scale Metabolic Modeling Elucidates the Role of 
# Proliferative Adaptation in Causing the Warburg Effect.")
# https://www.ebi.ac.uk/biomodels/MODEL1105100000
# https://www.ebi.ac.uk/compneur-srv/biomodels-main/MODEL1105100000

module MODEL1105100000
    import ..Chemostat_Rath2017: PROJ_ROOT, PROCESSED_DATA_DIR, RathData
    import Chemostat
    Ch = Chemostat
    import CSV
    import Serialization: deserialize
    import DataFrames: DataFrame
    

    include("meta.jl")
    include("dir_and_files.jl")
    include("load_data.jl")
    include("base_intake_info.jl")
    include("add_a1at_syntesis.jl")
end