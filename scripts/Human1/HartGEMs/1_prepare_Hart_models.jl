using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using ProgressMeter
    using Base.Threads

    import Chemostat
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
# I will use the HumanGEM base model configuration as a template
# It assumes that HumanGEM is a superset of the contextualized Hart models

## ---------------------------------------------------------------------
# base models
let
    modelid = "base"
    for tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test

        hart_model = AG.load_model(:raw, tissue)
        hart_model = MetNets.force_dims(hart_model)

        hart_model = ChR.prepare_metnet(AG, hart_model);
        
        println()
        @info("Saving", modelid, tissue)
        AG.save_model(modelid, tissue, hart_model)

        println()
        @info("Done!", tissue)

        println()
    end # for tissue
end

## ---------------------------------------------------------------------
# scaled models
let
    modelid = "scaled"
    lift_bound = 1e8
    scale_factor = 50.0
    sglob(AG; lift_bound, scale_factor)

    for tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test
        
        println("-"^50)
        @info("Scaling", tissue, 
            scale_factor, 
            lift_bound
        )
        
        base_model = AG.load_model(:base, tissue; uncompress = true)
        @info("Processing", tissue, 
            base_model = size(base_model),
            nzrange = Chemostat.Utils.nzabs_range(base_model.S)
        )
        scaled_model = Ch.Utils.well_scaled_model(base_model, scale_factor; lift_bound)
        @info("Done", tissue, 
            scaled_model = size(scaled_model), 
            nzrange = Chemostat.Utils.nzabs_range(scaled_model.S)
        )

        println()
        @info("Comparing with experiments")
        ChR.compare_with_experimets(scaled_model)
        
        println()
        @info("Saving", tissue)
        AG.save_model(modelid, tissue, scaled_model)
        println()
    
    end # for tissue
end

## ---------------------------------------------------------------------
# Fva models
let
    for baseid in ["base", "scaled"], tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test

        base_model = AG.load_model(baseid, tissue; uncompress = true)
        modelid = "fva_$(baseid)"
        
        println("-"^50)
        @info("FVA processing", modelid, tissue, 
            base_model = size(base_model),
        )

        # This run in parallel
        fva_model = Chemostat.LP.fva_preprocess(
            base_model; 
            check_obj = HG.HUMAN_BIOMASS_IDER
        )

        println()
        @info("Done!!", tissue, modelid,
            fva_model = size(fva_model)
        )

        println()
        @info("Comparing with experiments")
        ChR.compare_with_experimets(fva_model)
        
        println()
        @info("Saving", modelid, tissue)
        AG.save_model(modelid, tissue, fva_model)

        println()
        @info("Done!", modelid, tissue)

        println()

    end
end

