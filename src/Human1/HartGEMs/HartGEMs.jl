module HartGEMs

    import ...Chemostat_Rath2017
    import ..Human1
    import MAT
    import Chemostat
    import Chemostat.MetNets
    using ProjAssistant
    @gen_sub_proj

    include("const.jl")
    include("metnet.jl")

    function __init__()
        @create_proj_dirs
    end
end