# make 'fba' tests to all given lbs
function _is_essential!(lp_model::JuMP.Model, ridx::Int;
        test_lbs::Vector, 
        refidx::Integer, 
        refval::Float64, 
        var_th = 0.1,
        fba::Function, 
    )

    for lb_ in test_lbs
        MetLP.lb!(lp_model, ridx, lb_)
        fbaout = fba(lp_model)
        objval = MetLP.av(fbaout, refidx)
        abs(objval) > var_th * abs(refval) && return false
    end
    return true
end

function find_glc_min_glc_essential_medium!(model::MetNets.MetNet, biomidx::Int, exglcidx::Int;
        clear_caches = false,
        tries = 10_000,
        var_th = 0.1,
        solver = Clp.Optimizer,
        fba::Function = (model_) -> MetLP.fba(model_, biomidx)
    )

    # maps
    met_readable = HG.load_met_readable_ids()
    met_map = HG.load_exch_met_map()

    # find open
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

    # biom0
    fbaout0 = fba(model)
    biom0 = MetNets.av(model, fbaout0, biomidx)
    @show biom0

    # find assentials
    @info("Finding assentials")
    cid = (:ESSENTIALS, exchis)
    clear_caches && delcache(cid) # Reset
    essentialis = lcache(AG, cid) do
        _essentialis = Int[]
        for exchi in exchis
            lb_ = MetNets.lb(model, exchi)
            MetNets.lb!(model, exchi, 0.0)
            
            fbaout = fba(model)
            biom = MetNets.av(model, fbaout, biomidx)
            (biom == 0.0) && push!(_essentialis, exchi)

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
    fbaout = fba(model)
    biom =  MetNets.av(model, fbaout, biomidx)
    @show biom

    @info("Test adding glucose")
    MetNets.lb!(model, exglcidx, -1000.0)
    fbaout = fba(model)
    biom =  MetNets.av(model, fbaout, biomidx)
    @show biom

    # Experimental maximal "GLC" intake
    MAX_GLC = -(Rd.qval(:GLC) .|> abs |> maximum )
    @show MAX_GLC

    # I will find the minimum medium with make glucose essential
    @threads for _ in collect(1:tries)
        thid = threadid()
        
        # local copy
        # model_th = MetNets.MetNet(model; 
        #     lb = copy(model.lb), 
        #     ub = copy(model.ub), 
        #     c = copy(model.c)
        # )

        # lp_model
        lp_model = MetLP.build_lp_model(model; solver)
        # MetLP.lb!(lp_model, exchi, lb_)

        @info("Finding a greedy minimum medium", thid)
        minimum_mediumis = Set(exchis)
        for it in 1:(2*length(exchis))
            
            # Open minimum_mediumis
            for exchi in exchis
                # MetNets.lb!(model_th, exchi, -0.0)
                MetLP.lb!(lp_model, exchi, -0.0)
            end
            for exchi in minimum_mediumis
                # MetNets.lb!(model_th, exchi, -1000.0)
                MetLP.lb!(lp_model, exchi, -1000.0)
            end
            
            # pick one (not gluc)
            rexchi = rand(setdiff(minimum_mediumis, essentialis))
            exch = model.rxns[rexchi]
            rmet = get(met_map, exch, exch)
            rmet = get(met_readable, rmet, rmet)
            rmet == "glucose[s]" && continue # protect glucose
            
            # If rexch is not essential 
            if !_is_essential!(lp_model, rexchi;
                test_lbs = [0.0], # only test full closing
                refidx = biomidx, 
                refval = biom0, 
                var_th, fba
            )

                @info("Met NOT essential", it, length(minimum_mediumis), rmet, biom0, thid)

                # remove
                delete!(minimum_mediumis, rexchi)

                # Test glc essensiality
                if _is_essential!(lp_model, exglcidx;
                        test_lbs = [0.5, 0.3, 0.0] .* MAX_GLC,
                        refidx = biomidx, 
                        refval = biom0, 
                        var_th, fba
                    )
                    
                    @info("GLC essential", it, length(minimum_mediumis), biom0, biom, thid)
                    minimum_medium_ids = model.rxns[collect(minimum_mediumis)]
                    @show minimum_medium_ids

                    # save
                    medium_hash = hash((:MIN_MEDIUM, minimum_medium_ids))
                    sdat(AG, minimum_medium_ids, 
                        "minimum_medium", "glc_essential", (;hash = medium_hash), ".jls";
                        verbose = true
                    )
                    break
                else
                    # reset exglcidx
                    MetLP.lb!(lp_model, exglcidx, -1000.0)
                end
            else
                # reset rexchi
                MetLP.lb!(lp_model, rexchi, -1000.0)
                continue
            end
        end
        
    end # threads

    return
end
