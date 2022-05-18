using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

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

    using Plots

    using Statistics
end

## ---------------------------------------------------------------------
# DESCRIPTION
# This converge MaxEnt assuming the culture is glucose limited

## ---------------------------------------------------------------------
let
    mets_map = HG.load_mets_map()
    for rath_met in ["GLC", "GLN", "GAL"]
        model_met = mets_map[rath_met]
        @show rath_met model_met
    end
end