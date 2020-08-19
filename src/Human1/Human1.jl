module Human1
    import Chemostat.Utils: rxnindex, metindex, lb!, lb, ub!, S, S!, b,
                            ub, del_met, rxn_mets, rxn_reacts,
                            del_rxn, rxn_str, bounds, bounds!,
                            MetNet, FBAout, met_rxns, isrev
    import Chemostat.LP: fba
    import Chemostat.SteadyState: apply_bound!

    include("HumanGEM/HumanGEM.jl")
    include("tINIT_GEMs/tINIT_GEMs.jl")
    include("Rep_Human1/Rep_Human1.jl")
    include("ecGEMs/ecGEMs.jl")
    include("Utils/Utils.jl")
end