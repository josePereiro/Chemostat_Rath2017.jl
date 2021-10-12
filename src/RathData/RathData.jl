
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
    
    include("const.jl")
    include("raw_dat.jl")
    include("proc_dat.jl")
    include("a1at_aa_rel_ab.jl")
    include("data_interface.jl")
    
    function __init__()
        @create_proj_dirs
        
        _load_a1at_aa_rel_ab()
        _load_rath_bundle()
        _define_interface()
    end

end