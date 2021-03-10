function build_ecModel(srcmodel, ec_refdata_vec::Vector; add_protless = true)
    mM, mN = size(srcmodel) .* (3,5) # heuristic
    
    merge_model = ChU.expanded_model(srcmodel, mM, mN)
    
    rxn_ecmaps = [ec_refdata[ECMAP_KEY] for ec_refdata in ec_refdata_vec]
    allowed_protless = union([ec_refdata[PROTLESS_KEY] for ec_refdata in ec_refdata_vec]...)
    kin_stoi = maximum(median.([ec_refdata[PROT_STOIS_KEY].kin_stois for ec_refdata in ec_refdata_vec]))
    draw_stoi = maximum(median.([ec_refdata[PROT_STOIS_KEY].draw_stois for ec_refdata in ec_refdata_vec]))


    # Add extras
    for rxn_ecmap in rxn_ecmaps
        for rxndata in rxn_ecmap[EXTRAS_KEY]
            merge_rxndata!(merge_model, rxndata)
        end
    end
    
    srcrxns_count = length(srcmodel.rxns)
    prog = Progress(srcrxns_count; desc = "Merge ecmaps data   ")
    
    # Add rxn_ecmaps data
    for (srcrxni, srcrxn) in srcmodel.rxns |> enumerate
        for rxn_ecmap in rxn_ecmaps
            if haskey(rxn_ecmap, srcrxn)
                src_rxndata = RxnData(merge_model, srcrxni)
                clear_rxndata!(merge_model, src_rxndata)
                for rxndata in rxn_ecmap[srcrxn]
                    merge_rxndata!(merge_model, rxndata)
                end
                break
            end
        end
        update!(prog, srcrxni)
    end
    
    # Add protless
    if add_protless
        println("add_protless = true") # Test
        intersect!(allowed_protless, merge_model.rxns)
        merge_protless!(merge_model, allowed_protless, kin_stoi, draw_stoi);
    end

    finish!(prog)
    return ChU.compacted_model(merge_model)
end

function merge_rxndata!(model, rxndata)
    # check if already exist
    (rxndata.rxn in model.rxns) && return model
    
    # get empty spot
    nrxni = ChU.findempty(model, :rxns)
    
    # add data
    model.rxns[nrxni] = rxndata.rxn
    model.lb[nrxni] = rxndata.lb
    model.ub[nrxni] = rxndata.ub
    model.subSystems[nrxni] = rxndata.subSys
    
    for (i, met) in rxndata.mets |> enumerate
        # check if already exist
        if (met in model.mets)
            nmeti = ChU.metindex(model, met)
        else
            # get empty spot
            nmeti = ChU.findempty(model, :mets)
        end
        
        # add data
        model.mets[nmeti] = rxndata.mets[i]
        model.S[nmeti, nrxni] = rxndata.stoi[i]
        model.b[nmeti] = rxndata.b[i]
    end
    
    return model
end

function clear_rxndata!(model, rxndata)    
    # find reaction
    nrxni = ChU.rxnindex(model, rxndata.rxn)
    
    # clear data
    model.rxns[nrxni] = EMPTY_SPOT
    model.lb[nrxni] = 0.0
    model.ub[nrxni] = 0.0
    model.subSystems[nrxni] = "NA"
    
    for met in rxndata.mets
        nmeti = ChU.metindex(model, met)

        # clear stoi
        model.S[nmeti, nrxni] = 0.0

        # Do not delete if it is used in other reactions
        rxns = ChU.met_rxns(model, nmeti)
        length(rxns) > 0 && continue 
        model.mets[nmeti] = EMPTY_SPOT
        model.b[nmeti] = 0.0
    end
    
    return model
end

function print_ec_stats(ec_model)
    M, N = size(ec_model)
    # protless reactions
    protless_count = N
    draw_count = 0
    for (rxni, rxn) in ec_model.rxns |> enumerate
        metis = ChU.rxn_mets(ec_model, rxni)
        any(startswith(ec_model.mets[meti], PROT_PREFFIX) 
            for meti in metis) && (protless_count -= 1)
        startswith(rxn, DRAW_PREFFIX) && (draw_count += 1)
    end
    println("protless rxns  : ", protless_count, "/", N, " [", (protless_count * 100) ÷ N, "%]")
    println("draw rxns      : ", draw_count, "/", N, " [", (draw_count * 100) ÷ N, "%]")
    
    
    prot_met_count = 0
    pmet_count = 0
    for met in ec_model.mets
        startswith(met, PROT_PREFFIX) && (prot_met_count += 1)
        startswith(met, PMET_PREFFIX) && (pmet_count += 1)
    end
    println("prot mets      : ", prot_met_count, "/", M, " [", (prot_met_count * 100) ÷ M, "%]")
    println("pmets          : ", pmet_count, "/", M, " [", (pmet_count * 100) ÷ M, "%]")
    
end