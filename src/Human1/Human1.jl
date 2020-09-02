module Human1
    
    import Chemostat.Utils: rxnindex, metindex, lb!, lb, ub!, S, S!, b,
                            ub, rxn_mets, rxn_reacts, isfwd,
                            rxn_str, bounds, bounds!, is_exchange,
                            MetNet, FBAout, met_rxns, isrev, exchanges, 
                            summary
    import Chemostat.LP: fba
    import Chemostat.SteadyState: apply_bound!
    using ProgressMeter
    import Distributions: median, mean
    import ..Chemostat_Rath2017: RathData

    include("HumanGEM/HumanGEM.jl")
    include("tINIT_GEMs/tINIT_GEMs.jl")
    include("Rep_Human1/Rep_Human1.jl")
    include("ecGEMs/ecGEMs.jl")
    include("Utils/Utils.jl")

    function __init__()
        # Utils
        _ider_maps()
    end
end