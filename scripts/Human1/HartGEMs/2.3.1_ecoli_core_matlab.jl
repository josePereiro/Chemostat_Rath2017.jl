using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using Plots
    import GR
    isinteractive() && GR.inline("png")
    using Statistics
    using ProgressMeter
    using Base.Threads

    import Chemostat
    import Chemostat.MetNets
    import Chemostat.MetLP
    import Chemostat.MetEP
    const Ch = Chemostat

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
    const AG = H1.HartGEMs

    using MATLAB

end

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
# setup CobraToolbox
let
    mat"""
        addpath('/Users/Pereiro/University/Research/Metabolism/MaxEntEP2020/WIP/Chemostat_Human1/deps/cobratoolbox');
        initCobraToolbox(false)
    """
    # mat"disp('Hii')"
    # MetNets.ECOLI_MODEL_MAT_FILE
end

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
include("2.0.3_glc_gradual_limiting_check.jl") 
let
    # prepare model
    model = MetNets.ecoli_core_model()

    biomidx = MetNets.rxnindex(model, MetNets.ECOLI_MODEL_BIOMASS_IDER)
    exglcidx = MetNets.rxnindex(model, "EX_glc(e)")

    for solver_mod in [Clp, Cbc, GLPK]
        solver = solver_mod.Optimizer
        head_hits = ["minimum_medium", "glc_essential"]
        for name in collect(readdir(datdir(AG)))
            heads, params, _ = parse_dfname(name)
            heads != head_hits && continue
            
            println("\n", "-"^60, "\n")
            @show name
            file = joinpath(datdir(AG), name)
            medium_ids = ldat(AG, file)

            p = plot_glc_gradual_limiting_check(
                model, medium_ids, biomidx, exglcidx;
                exglclb0 = -1000.0,
                fba = (model_) -> MetLP.fba(model_, biomidx),
                fs = [0.0; 10.0 .^ (-3:0.05:0.0)],
                # fs = [0.0, 0.5, 1.0], 
                solver
            )

            sfig(AG, p, heads..., params, (;solver = string(solver_mod)), ".png"; verbose = true)
        end
    end
end