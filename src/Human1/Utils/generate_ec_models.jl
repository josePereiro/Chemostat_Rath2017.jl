function build_ecModel(srcmodel, ec_refdatas::Vector)
    mM, mN = size(srcmodel) .* (3,5) # heuristic
    
    merge_model = expanded_model(srcmodel, mM, mN)
    
    rxn_ecmaps = [ec_refdata[ECMAP_KEY] for ec_refdata in ec_refdatas]
    allowed_protless = union([ec_refdata[PROTLESS_KEY] for ec_refdata in ec_refdatas]...)
    kin_stoi = maximum(median.([ec_refdata[PROT_STOIS_KEY].kin_stois for ec_refdata in ec_refdatas]))
    draw_stoi = maximum(median.([ec_refdata[PROT_STOIS_KEY].draw_stois for ec_refdata in ec_refdatas]))


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
    intersect!(allowed_protless, merge_model.rxns)
    merge_protless!(merge_model, allowed_protless, kin_stoi, draw_stoi);

    finish!(prog)
    return compated_model(merge_model)
#     return (merge_model, allowed_protless, kin_stoi, draw_stoi)
end


# This just prepare the model to hold more elements by making 
# all the numerical fields larger
function expanded_model(model, expdM::Int, expdN::Int)
    M, N = size(model)
    @assert all((expdM, expdN) .> (M, N))
    
    S_ = zeros(expdM, expdN)
    S_[1:M, 1:N] .= model.S
    b_ = zeros(expdM)
    b_[1:M] .= model.b
    lb_ = zeros(expdN)
    lb_[1:N] .= model.lb
    ub_ = zeros(expdN)
    ub_[1:N] .= model.ub
    subSystems_ = Vector{Any}(fill("NA", expdN))
    subSystems_[1:N] .= model.subSystems
    rxns_ = fill(EMPTY_SPOT, expdN)
    rxns_[1:N] .= model.rxns
    mets_ = fill(EMPTY_SPOT, expdM)
    mets_[1:M] .= model.mets
    
    return MetNet(S_, b_, lb_, ub_, rxns_, mets_; subSystems = subSystems_)
end

function free_spot(model, col)
    spot = findfirst(isequal(EMPTY_SPOT), getfield(model, col))
    isnothing(spot) && error("Not $col empty spot!!!")
    return spot
end

function compated_model(model)

    empty_mets = findall(model.mets .== EMPTY_SPOT)
    met_idxs = trues(size(model, 1))
    met_idxs[empty_mets] .= false
    M = count(met_idxs)

    empty_rxns = findall(model.rxns .== EMPTY_SPOT)
    rxn_idxs = trues(size(model, 2))
    rxn_idxs[empty_rxns] .= false
    N = count(rxn_idxs)
    
    S_ = model.S[met_idxs, rxn_idxs]
    b_ = model.b[met_idxs]
    lb_ = model.lb[rxn_idxs]
    ub_ = model.ub[rxn_idxs]
    subSystems_ = model.subSystems[rxn_idxs]
    rxns_ = model.rxns[rxn_idxs]
    mets_ = model.mets[met_idxs]
    
    return MetNet(S_, b_, lb_, ub_, rxns_, mets_; subSystems = subSystems_)
    
end


function merge_rxndata!(model, rxndata)
    # check if already exist
    (rxndata.rxn in model.rxns) && return model
    
    # get empty spot
    nrxni = free_spot(model, :rxns)
    
    # add data
    model.rxns[nrxni] = rxndata.rxn
    model.lb[nrxni] = rxndata.lb
    model.ub[nrxni] = rxndata.ub
    model.subSystems[nrxni] = rxndata.subSys
    
    for (i, met) in rxndata.mets |> enumerate
        # check if already exist
        if (met in model.mets)
            nmeti = metindex(model, met)
        else
            # get empty spot
            nmeti = free_spot(model, :mets)
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
    nrxni = free_spot(model, :rxns)
    
    # clear data
    model.rxns[nrxni] = EMPTY_SPOT
    model.lb[nrxni] = 0.0
    model.ub[nrxni] = 0.0
    model.subSystems[nrxni] = "NA"
    
    for met in rxndata.mets
        nmeti = free_spot(model, :mets)

        # clear stoi
        model.S[nmeti, nrxni] = 0.0

        # Do not delete if it is used in other reactions
        rxns = met_rxns(model, nmeti)
        length(rxns) > 0 && continue 
        model.mets[nmeti] = EMPTY_SPOT
        model.b[nmeti] = 0.0
    end
    
    return model
end
