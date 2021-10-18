module Human1

    using ProjAssistant
    @gen_sub_proj

    include("HumanGEM/HumanGEM.jl")
    include("HartGEMs/HartGEMs.jl")

    function __init__()
        @create_proj_dirs
    end

end