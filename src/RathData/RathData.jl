
# This module interface the data retrived from Rath, Alexander. 
# “Characterisation of Cell Growth, Metabolism and 
# Recombinant Protein Production during Transient and 
# Steady State Conditions for the Human Cell Line AGE1.HN-AAT,” 
# 2017. https://pure.mpg.de/pubman/item/item_2508673_4.

module RathData

    import CSV

    using DataFrames
    using Serialization
    using ProjAssistant
    @gen_sub_proj
    
    include("iders.jl")
    # include("dir_and_files.jl")
    # include("data_interface.jl")
    # include("a1at_aa_rel_ab.jl")
    
    function __init__()
        @create_proj_dirs
        # _load_rath_bundle()
        # _define_interface()
    end

end