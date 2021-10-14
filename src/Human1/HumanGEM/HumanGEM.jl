module HumanGEM

    import ...Chemostat_Rath2017
    import ..Human1
    import MAT
    import Chemostat
    using ProjAssistant
    @gen_sub_proj

    include("const.jl")
    include("mets_map.jl")
    include("Hams_medium.jl")
    include("metnet.jl")
    include("Niklas_biomass.jl")
    include("base_intake_info.jl")

    function __init__()
        @create_proj_dirs
    end
end
