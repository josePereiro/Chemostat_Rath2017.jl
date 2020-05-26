
# This module interface the data retrived from from Rath, Alexander. 
# “Characterisation of Cell Growth, Metabolism and 
# Recombinant Protein Production during Transient and 
# Steady State Conditions for the Human Cell Line AGE1.HN-AAT,” 
# 2017. https://pure.mpg.de/pubman/item/item_2508673_4.

module RathData

    import ..Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    using DataFrames
    using Serialization
    import CSV

    include("iders.jl")
    include("dir_and_files.jl")
    include("data_interface.jl")
    include("a1at_aa_rel_ab.jl")
end