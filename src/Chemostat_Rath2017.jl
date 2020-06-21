"""
    This package contains both, the script that will produce results
    and the common tools used for produced.
    The module itself only contains the common tools
"""
module Chemostat_Rath2017

    include("Utils/Utils.jl")
    include("RathData/RathData.jl")
    include("MODEL1105100000/MODEL1105100000.jl")
    include("Human1/Human1.jl")

end # module
