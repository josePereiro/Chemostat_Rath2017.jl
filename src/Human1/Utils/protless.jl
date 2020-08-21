function collect_protless(model::MetNet, idxs = eachindex(model.rxns))
    filter(idxs) do rxni
        metis = rxn_mets(model, rxni)
        metids = model.mets[metis]
        !any(startswith(met, PROT_PREFFIX) for met in metids)
    end
end

function prot_stois(ref_model)
    prot_kin_stois = []
    prot_draw_stois = []

    prot_pool_idx = metindex(ref_model, PROT_POOL)
    foreach(eachindex(ref_model.mets)) do meti
        met = ref_model.mets[meti]
        !startswith(met, PROT_PREFFIX) && return
        meti == prot_pool_idx && return

        rxnis = met_rxns(ref_model, meti)
        for rxni in rxnis
            rxn = ref_model.rxns[rxni]
            if startswith(rxn, DRAW_PREFFIX)
                s = S(ref_model, prot_pool_idx, rxni)
                push!(prot_draw_stois, s) 
            else
                s = S(ref_model, meti, rxni)
                push!(prot_kin_stois, s)
            end
        end
    end

    return (kin_stois = prot_kin_stois, draw_stois = prot_draw_stois)
end

function merge_protless!(model, allowed_protless, prot_kin_stoi, prot_draw_stoi)
    
    @assert prot_kin_stoi < 0.0
    @assert prot_draw_stoi < 0.0
    
    all_protless = model.rxns[collect_protless(model)]
    disallowed_protless = filter(setdiff(all_protless, allowed_protless)) do rxn
        rxn != EMPTY_SPOT && !is_exchange(model, rxn)
    end |> unique!
    
    disallowed_protless_count = length(disallowed_protless)
    prog = Progress(disallowed_protless_count; desc = "Merge protless  ")
    for (i, fwd_rxn) in disallowed_protless |> enumerate
        
        fwd_rxni = rxnindex(model, fwd_rxn) # forward defined
        rev = isrev(model, fwd_rxni)
        
        ## forward reaction
        # A + (-prot_kin_stoi) prot_met -> B
        # add prot_met
        prot_met = "prot_" * fwd_rxn # This must make it unique
        prot_meti = free_spot(model, :mets)
        model.S[prot_meti, fwd_rxni] = prot_kin_stoi
        model.mets[prot_meti] = prot_met
        model.b[prot_meti] = 0.0
        
        # Adding draw rxn
        # (-prot_draw_stoi) prot_pool ==> (1.0) prot_met
        draw_rxni = free_spot(model, :rxns)
        model.rxns[draw_rxni] = DRAW_PREFFIX * prot_met
        model.lb[draw_rxni], model.ub[draw_rxni] = (0.0, 1000.0)
        model.S[prot_meti, draw_rxni] = 1.0
        prot_pooli = metindex(model, Human1.PROT_POOL)
        model.S[prot_pooli, draw_rxni] = prot_draw_stoi
        model.subSystems[draw_rxni] = model.subSystems[fwd_rxni]
        
        if rev
            ## backward reaction
            # B + (-prot_kin_stoi) prot_met -> A
            bkwd_rxni = free_spot(model, :rxns)
            model.S[:, bkwd_rxni] .= -model.S[:, fwd_rxni]
            model.S[prot_meti, bkwd_rxni] = prot_kin_stoi
            model.lb[bkwd_rxni] = 0.0
            model.ub[bkwd_rxni] = abs(model.lb[fwd_rxni])
            model.rxns[bkwd_rxni] = fwd_rxn * REV_SUFFIX
            model.subSystems[bkwd_rxni] = model.subSystems[fwd_rxni]
            
            # make forward irrev
            model.lb[fwd_rxni] = 0.0 
        end
        
        # progress
        update!(prog, i)
    end
    # progress
    finish!(prog)
    
    return model
    
end