
function fill_protless(model::MetNet, protless_rxns::Vector, 
    prot_kin_stoi, prot_draw_stoi)
    @assert prot_kin_stoi < 0.0
    @assert prot_draw_stoi < 0.0

    # Computing final dimentions
    nM, nN = size(model)
    println("Current dimention: ", (nM, nN))
    for rxni in protless_rxns
        rev = isrev(model, rxni)
        # One new reaction ('_REV') for each rev
        rev && (nN += 1)
        # One new prot_ met for each protless reaction
        nM += 1
        # One new 'draw_' reaction for each 'preot_' met
        nN += 1

    end
    println("Final dimention: ", (nM, nN))

    M, N = size(model)
    S_ = zeros(nM, nN);
    S_[1:M, 1:N] .= model.S
    b_ = zeros(nM);
    b_[1:M] .= model.b
    lb_ = zeros(nN);
    lb_[1:N] .= model.lb
    ub_ = zeros(nN);
    ub_[1:N] .= model.ub
    rxns_ = Vector{String}(undef, nN);
    rxns_[1:N] .= model.rxns
    mets_ = Vector{String}(undef, nM);
    mets_[1:M] .= model.mets
    subSystems_ = Vector(undef, nN);
    subSystems_[1:N] .= model.subSystems;


    nmeti = size(model, 1) + 1
    nrxni = size(model, 2) + 1
    # This need to be stimated from template
    prot_pool_idx = metindex(model, PROT_POOL)

    for rxni in protless_rxns
        rev = isrev(model, rxni)
        rxn = model.rxns[rxni]

        ## forward defined
        # add prot_met
        S_[nmeti, rxni] = prot_kin_stoi
        mets_[nmeti] = "prot_" * rxn
        lb_[rxni] = rev ? 0.0 : lb_[rxni]
        b_[nmeti] = 0.0

        # Adding draw rxn
        rxns_[nrxni] = "draw_" * mets_[nmeti]
        lb_[nrxni], ub_[nrxni] = (0.0, 1000.0)
        # (-prot_draw_stoi) prot_pool ==> (1.0) prot_met
        S_[nmeti, nrxni] = 1.0
        S_[prot_pool_idx, nrxni] = prot_draw_stoi
        subSystems_[nrxni] = subSystems_[rxni]
        nmeti += 1  
        nrxni += 1

        ## backward defined
        # prot_met was added in the foward reaction
        if rev
            S_[:, nrxni] .= -S_[:, rxni]
            lb_[nrxni] = 0.0
            ub_[nrxni] = abs(lb_[rxni])
            rxns_[nrxni] = rxns_[rxni] * REV_SUFFIX
            subSystems_[nrxni] = subSystems_[rxni]
            nrxni += 1
        end

    end

    return MetNet(S_, b_, lb_, ub_, rxns_, mets_; subSystems = subSystems_);
end

