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
# # Load model
# let
#     fn = "/Users/Pereiro/University/Research/Metabolism/MaxEntEP2020/WIP/Chemostat_Rath2017/data/raw/Human1/Human1_Publication_Data_Scripts/ec_GEMs/models/HOP62/ecModel_batch.mat"
#     ec_model = MetNets.read_mat(fn)
#     MetNets.clampfields!(ec_model, [:lb, :ub]; abs_max = 999999, zeroth = 1e-8)
#     objidx = MetNets.rxnindex(ec_model, HG.HUMAN_BIOMASS_IDER)
#     ec_model = ChR.prepare_metnet(AG, ec_model; inf_medium = true);
#     ec_model = MetNets.force_dims(ec_model)
#     # fbaout = MetLP.fba!(ec_model, objidx)
#     sdat(AG, MetNets.compressed_model(ec_model), "ec_model", ".jls")
# end 

## ---------------------------------------------------------------------
let

    met_readable = HG.load_met_readable_ids()
    met_map = HG.load_exch_met_map()

    @info("Loading model")
    # model = MetNets.uncompressed_model(ldat(AG, "ec_model", ".jls"))
    tissue = "GBM"
    modelid = "base"
    model = AG.load_model(;modelid, tissue, uncompress = true)
    objidx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER)
    exglcidx = MetNets.rxnindex(model, HG.HUMAN_GLC_EX_IDER)

    @info("Finding open exchanges")
    exchis = filter(MetNets.exchanges(model)) do exchi
        MetNets.lb(model, exchi) < 0
    end
    @show length(exchis)
    
    # Open all
    @info("Fully opening all")
    for exchi in exchis
        MetNets.lb!(model, exchi, -1000)
    end
    fbaout = MetLP.fba!(model, objidx)
    biom0 = MetLP.objval(fbaout)
    @show biom0

    # find assentials
    @info("Finding assentials")
    cid = (:ESSENTIALS, exchis)
    # delcache(cid) # Reset
    essentialis = lcache(AG, cid) do
        _essentialis = Int[]
        for exchi in exchis
            lb_ = MetNets.lb(model, exchi)
            MetNets.lb!(model, exchi, 0.0)
            
            fbaout = MetLP.fba!(model, objidx)
            biom = MetLP.objval(fbaout)
            (biom == 0.0) && push!(_essentialis, exchi)

            exch = model.rxns[exchi]
            met = get(met_map, exch, exch)
            met = get(met_readable, met, met)

            @info(met, biom)

            MetNets.lb!(model, exchi, lb_)
        end
        return _essentialis
    end
    @show length(essentialis)

    @info("non_essentialis")
    non_essentialis = setdiff(exchis, essentialis)
    @show length(non_essentialis)

    @info("close non_essentials")
    for exchi in exchis
        MetNets.lb!(model, exchi, -1000.0)
    end
    for exchi in non_essentialis
        MetNets.lb!(model, exchi, 0.0)
    end
    
    @info("Test only essentials")
    fbaout = MetLP.fba!(model, objidx)
    biom = MetLP.objval(fbaout)
    @show biom

    @info("Test adding glucose")
    MetNets.lb!(model, exglcidx, -1000.0)
    fbaout = MetLP.fba!(model, objidx)
    biom = MetLP.objval(fbaout)
    @show biom

    # I will find the minimum medium with make glucose essential
    @threads for _ in collect(1:10000)
        thid = threadid()
        
        # local copy
        model_th = MetNets.MetNet(model; 
            lb = copy(model.lb), 
            ub = copy(model.ub), 
            c = copy(model.c)
        )

        @info("Finding a greedy minimum medium", thid)
        minimum_mediumis = Set(exchis)
        for it in 1:(3*length(exchis))
            
            # Open minimum_mediumis
            for exchi in exchis
                MetNets.lb!(model_th, exchi, -0.0)
            end
            for exchi in minimum_mediumis
                MetNets.lb!(model_th, exchi, -1000.0)
            end
            MetNets.lb!(model_th, exglcidx, -1000.0)
            
            # close one
            rexchi = rand(setdiff(minimum_mediumis, essentialis))
            exch = model_th.rxns[rexchi]
            rmet = get(met_map, exch, exch)
            rmet = get(met_readable, rmet, rmet)
            rmet == "glucose[s]" && continue # protect glucose
            MetNets.lb!(model_th, rexchi, 0.0)
            
            fbaout = MetLP.fba!(model_th, objidx)
            biom = MetLP.objval(fbaout)
            
            # If rexch is not essential 
            if abs(biom) > 0.9 * abs(biom0)

                @info("Met NOT essential", it, length(minimum_mediumis), rmet, biom0, biom, thid)

                # remove
                delete!(minimum_mediumis, rexchi)

                # Test glc essensiality
                MetNets.lb!(model_th, exglcidx, 0.0)
                fbaout = MetLP.fba!(model_th, objidx)
                biom = MetLP.objval(fbaout)
                # If glc is essential
                if abs(biom) < 0.1 * abs(biom0)

                    @info("GLC essential", it, length(minimum_mediumis), biom0, biom, thid)
                    @show minimum_mediumis

                    # save
                    medium_hash = hash((:MIN_MEDIUM, minimum_mediumis))
                    sdat(AG, minimum_mediumis, 
                        "minimum_medium", (;medium_hash), ".jls";
                        verbose = true
                    )
                else
                    @info("GLC NOT essential", it, length(minimum_mediumis), biom0, biom, thid)
                end
            else
                # If rexch is essential (open)
                @info("Met essential", it, length(minimum_mediumis), rmet, biom0, biom, thid)
                MetNets.lb!(model_th, rexchi, -1000.0)
                continue
            end
            println()
        end
        
    end # threads


end

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
# ## ---------------------------------------------------------------------
# # histogram
# let
#     dir = datdir(AG)
#     lens = Int[]
#     medium_files = filter(readdir(dir; join = true)) do fn
#         contains(fn, "minimum_medium")
#     end
#     @show length(medium_files)
#     for fn in medium_files
#         push!(lens, length(ldat(fn)))
#     end
#     histogram(lens; bins = 20)
# end


# ## ---------------------------------------------------------------------
# let
#     dir = datdir(AG)
#     lens = Int[]
#     medium_files = filter(readdir(dir; join = true)) do fn
#         contains(fn, "minimum_medium")
#     end

#     met_readable = HG.load_met_readable_ids()
#     met_map = HG.load_exch_met_map()

#     @info("Loading model")
#     # model = MetNets.uncompressed_model(ldat(AG, "ec_model", ".jls"))
#     tissue = "GBM"
#     modelid = "base"
#     model = AG.load_model(;modelid, tissue, uncompress = true)
#     objidx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER)

#     @info("Finding open exchanges")
#     exchis = filter(MetNets.exchanges(model)) do exchi
#         MetNets.lb(model, exchi) < 0
#     end
#     @show length(exchis)

#     # test
#     fbaout = MetLP.fba!(model, objidx)
#     biom0 = MetLP.objval(fbaout)
#     @show biom0

#     for medium_file in medium_files
        
#         println()
#         minimum_mediumis = ldat(medium_file)
#         @show minimum_mediumis
#         @show length(minimum_mediumis)
        
#         println()
#         @show biom0

#         # close all
#         for exchi in exchis
#             MetNets.lb!(model, exchi, 0.0)
#         end
        
#         # open minimum
#         for exchi in minimum_mediumis
#             MetNets.lb!(model, exchi, -1000.0)
#         end

#         # test
#         fbaout = MetLP.fba!(model, objidx)
#         biom = MetLP.objval(fbaout)
#         @show biom
        
#         # close glcose
#         MetNets.lb!(model, HG.HUMAN_GLC_EX_IDER, 0.0)
        
#         # test
#         fbaout = MetLP.fba!(model, objidx)
#         biom = MetLP.objval(fbaout)
#         @show biom

#     end

# end

# ## ---------------------------------------------------------------------
# ## ---------------------------------------------------------------------
# # let
# #     met_readable = HG.load_met_readable_ids()
# #     met_map = HG.load_exch_met_map()
# #     results = ldat(AG, "glc_limited_study", ".jls")
    
# #     li = 5
# #     biomis = Float64[]
# #     mets = String[]
# #     p = plot()
# #     for (exch, dat) in results
# #         met = get(met_map, exch, exch)
# #         met = get(met_readable, met, met)
# #         println("met: ", met)
# #         println("bioms: ", join(dat, ", "))
# #         println()
# #         # idx = findfirst(reverse(dat)) do biom
# #         #     biom < maximum(dat) * 0.9
# #         # end
# #         # idx = isnothing(idx) ? 0 : idx
# #         if any(iszero.(dat))
# #             # push!(biomis, idx)
# #             # push!(mets, met)
# #             plot!(p, dat; label = met)
# #         end
# #     end
# #     p
# #     # rang_ = 1:10
# #     # sis = sortperm(biomis; rev = true)
# #     # bar(mets[sis], biomis[sis]; xrotation = 45)
    
# #     # sfig(AG, ps,
# #     #     @fileid, "glc_limited_study", ".png";
# #     #     layout = (3,3)
# #     # )
# # end