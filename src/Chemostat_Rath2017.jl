"""
    This package contains both, the script that will produce results
    and the common tools used for produced.
    The module itself only contains the common tools
"""
module Chemostat_Rath2017

    using ProjAssistant
    @gen_top_proj

    import Chemostat
    const Ch = Chemostat
    import Chemostat.MetNets
    import Chemostat.MetEP
    import Chemostat.MetLP

    
    include("RathData/RathData.jl")
    const Rd = RathData
    
    include("Human1/Human1.jl")
    const H1 = Human1
    const HG = H1.HumanGEM
    
    include("Utils/Utils.jl")
    
    function __init__()
        @create_proj_dirs
    end
    
end # module
