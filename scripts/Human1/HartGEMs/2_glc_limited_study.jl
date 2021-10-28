using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using Plots
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

end

## ---------------------------------------------------------------------
# DESCRIPTION
# The intake of glucose must be the limiting external flux of the network

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    modelid = "fva_base"
    model = AG.load_model(;modelid, tissue, uncompress = true)
    objidx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER)

    exchis = filter(MetNets.exchanges(model)) do exchi
        MetNets.lb(model, exchi) < 0
    end
    @show length(exchis)
    
    # Open all but glc
    for exchi in exchis
        MetNets.lb!(model, exchi, -1000)
    end
    MetNets.lb!(model, HG.HUMAN_GLC_EX_IDER, -0.5)

    results = Dict()
    @showprogress for exchi in exchis
        
        exch = model.rxns[exchi]
        dat = get!(results, exch, Float64[])
        
        MetNets.summary(model, model.rxns[exchi])
        println()
        
        lb_ = MetNets.lb(model, exchi)

        for f in 0.0:0.1:1.0
            MetNets.lb!(model, exchi, lb_ * f)
            fbaout = MetLP.fba!(model, objidx)
            biom = MetLP.objval(fbaout)
            push!(dat, biom)
        end

        MetNets.lb!(model, exchi, lb_)
    end

    sdat(AG, results, 
        "glc_limited_study", ".jls";
        verbose = true
    )

end;

## ---------------------------------------------------------------------
let
    met_readable = HG.load_met_readable_ids()
    met_map = HG.load_exch_met_map()
    results = ldat(AG, "glc_limited_study", ".jls")
    
    li = 5
    bioms = Float64[]
    mets = String[]
    for (exch, dat) in results
        met = get(met_map, exch, exch)
        met = get(met_readable, met, met)
        push!(bioms, dat[li])
        push!(mets, met)
    end
    
    rang_ = 1:10
    sis = sortperm(bioms)
    bar(mets[sis][rang_], bioms[sis][rang_]; xrotation = 45)
    
    # sfig(AG, ps,
    #     @fileid, "glc_limited_study", ".png";
    #     layout = (3,3)
    # )
end