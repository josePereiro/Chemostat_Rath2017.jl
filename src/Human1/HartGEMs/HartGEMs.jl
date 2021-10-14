module HartGEMs

    import ...Chemostat_Rath2017
    import ..Human1
    using ProjAssistant
    @gen_sub_proj

    include("metnet.jl")

    function __init__()
        @create_proj_dirs
    end
end