
# This module interface the data retrived from from Rath, Alexander. 
# “Characterisation of Cell Growth, Metabolism and 
# Recombinant Protein Production during Transient and 
# Steady State Conditions for the Human Cell Line AGE1.HN-AAT,” 
# 2017. https://pure.mpg.de/pubman/item/item_2508673_4.

module RathData

    import ..Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import DataFrames: DataFrame
    import CSV

    include("1_iders.jl")
    include("2_files.jl")
    include("3_convert_rath_data.jl")
    include("4_interface.jl")
end