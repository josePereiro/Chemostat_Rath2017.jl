using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using ProgressMeter
    using Base.Threads

    import Chemostat
    import Chemostat.MetNets
    import Chemostat.MetLP
    import Chemostat.MetEP

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
    const AG = H1.HartGEMs
end

## ---------------------------------------------------------------------
# DESCRIPTION
# It use the contextualized Hart GEMs and produce a set of other GEMs
# contextualized using Rath data

## ---------------------------------------------------------------------
# base models
let
    modelid = "base"
    for tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test

        AG.check_modelfile(;modelid, tissue) && continue

        raw_model = AG.load_model(;modelid = :raw, tissue, uncompress = true)
        raw_model = MetNets.force_dims(raw_model)

        hart_model = ChR.prepare_metnet(AG, raw_model; inf_medium = true);
        hart_model = MetNets.force_dims(hart_model)
        
        println()
        @info("Saving", modelid, tissue)
        AG.save_model(hart_model; modelid, tissue)

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
        
        @info("Scaling", 
            modelid, tissue, 
            scale_factor, 
            lift_bound
        )
        AG.check_modelfile(;modelid, tissue) && continue

        base_model = AG.load_model(;modelid = "base", tissue, uncompress = true)
        @info("Processing", tissue, 
            base_model = size(base_model),
            nzrange = MetNets.nzabs_range(base_model.S)
        )
        scaled_model = MetNets.well_scaled_model(base_model, scale_factor; lift_bound)
        scaled_model = MetNets.force_dims(scaled_model)
        @info("Done", tissue, 
            scaled_model = size(scaled_model), 
            nzrange = MetNets.nzabs_range(scaled_model.S)
        )

        println()
        @info("Comparing with experiments")
        ChR.compare_with_experimets(scaled_model)
        
        println()
        @info("Saving", tissue)
        AG.save_model(scaled_model; modelid, tissue)
        println()
    
    end # for tissue
end

## ---------------------------------------------------------------------
# Fva models
let
    for baseid in ["base", "scaled"], tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test

        base_model = AG.load_model(;modelid = baseid, tissue, uncompress = true)
        modelid = "fva_$(baseid)"
        
        println("-"^50)
        @info("FVA processing", modelid, tissue, 
            base_model = size(base_model),
        )
        AG.check_modelfile(;modelid, tissue) && continue

        # This run in parallel
        fva_model = MetLP.fva_preprocess(
            base_model; 
            check_obj = HG.HUMAN_BIOMASS_IDER
        ) |> MetNets.force_dims

        println()
        @info("Done!!", tissue, modelid,
            fva_model = size(fva_model)
        )

        println()
        @info("Comparing with experiments")
        ChR.compare_with_experimets(fva_model)
        
        println()
        @info("Saving", modelid, tissue)
        AG.save_model(fva_model; modelid, tissue)

        println()
        @info("Done!", modelid, tissue)

        println()

    end
end

## ---------------------------------------------------------------------
# stst models
let
    for stst in Rd.STSTS, baseid in ["scaled"], tissue in AG.TISSUES
        (tissue != "GBM") && continue # Test
        (stst != "A") && continue # Test

        modelid = "stst_$(baseid)"
        fva_model = AG.load_model(;modelid = "fva_$(baseid)", tissue, uncompress = true)
        
        println("-"^50)
        @info("STST processing", 
            modelid, tissue, stst,
            fva_model = size(fva_model),
        )
        AG.check_modelfile(;modelid, tissue, stst) && continue

        # stst contextualization
        intake_info = HG.stst_base_intake_info(stst)
        xi = Rd.val(:ξ, stst)
        Chemostat.apply_bound!(fva_model, xi, intake_info;
            emptyfirst = true, ignore_miss = true
        )

        # This run in parallel
        stst_model = MetLP.fva_preprocess(fva_model; 
            check_obj = HG.HUMAN_BIOMASS_IDER
        ) |> MetNets.force_dims

        println()
        @info("Done!!", tissue, modelid,
            stst_model = size(stst_model)
        )

        println()
        @info("Comparing with experiments")
        ChR.compare_with_experimets(stst_model)
        
        println()
        @info("Saving", modelid, tissue)
        AG.save_model(stst_model; modelid, tissue, stst)

        println()
        @info("Done!", modelid, tissue, stst)

        println()
        return

    end
end