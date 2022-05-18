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
# This is a numeric test of MaxEnt EP at beta 0

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    for modelid in ["fva_base", "fva_scaled"]

        model = AG.load_model(modelid, tissue; uncompress = true)
        model = MetNets.force_dims(model)
        objider = HG.HUMAN_BIOMASS_IDER
        M, N = size(model)
        objidx = MetNets.rxnindex(model, objider)

        # maxent
        save_frec = 5
        function oniter(it, epmodel)
            !iszero(rem(it, save_frec)) && return
            
            epout_ = MetEP.produce_epout(epmodel; drop_epfields = true)
            sdat(AG, epout_,
                "maxent_ep_beta0", "epout", (;modelid, tissue, it), ".jls";
                verbose = false
            )
            
        end
        
        try
            @info("MAxEnt", model = size(model))
            epmodel = MetEP.EPModel(model)
            model = nothing; GC.gc(); # free mem
            MetEP.converge_ep!(epmodel; oniter,
                maxiter = 2000, # Test 
                verbose = true
            )
        catch err
            @error("At MaxEnt ", err)
        end
    end

end