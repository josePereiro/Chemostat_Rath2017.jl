"""
    This package contains both, the script that will produce results
    and the common tools used for produced.
    The module itself only contains the common tools
"""
module Chemostat_Rath2017

    using ProjAssistant
    @gen_top_proj

    include("RathData/RathData.jl")

    function __init__()
        @create_proj_dirs
    end
    
end # module
