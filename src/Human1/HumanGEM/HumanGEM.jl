module HumanGEM

    import ...Chemostat_Rath2017
    import ..Human1
    import MAT
    import Chemostat
    using ProjAssistant
    @gen_sub_proj

    include("base_intake_info.jl")
    include("const.jl")
    include("exch_met_map.jl")
    include("Hams_medium.jl")
    include("met_readable_ids.jl")
    include("metnet.jl")
    include("mets_map.jl")
    include("Niklas_biomass.jl")

    function __init__()
        @create_proj_dirs
    end
end
