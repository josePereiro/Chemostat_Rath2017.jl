# this dict contains all the match for each src rxn id found in
# the ec models
function get_rxn_ecmap(src_model, ec_model; verbose = true)
    
    dict = Dict{String, Vector{RxnData}}()
    ec_rxns_set = Set(ec_model.rxns)
    
    srcN = size(src_model, 2)
    ecN = size(ec_model, 2)
    count = 0
    
    # Info
    prog = Progress(srcN; desc = "Build ec_map   ")
    showvalues() = [("src matched     ", string(count, " [", (count * 100) ÷ srcN, "%]")),
                ("left in ec_model", string(length(ec_rxns_set), " [", (length(ec_rxns_set) * 100) ÷ ecN, "%]"))]
    
    # extract src -> ec relation
    for (i, srcrxn) in src_model.rxns |> enumerate
        
        dict[srcrxn] = RxnData[]
        
        # Pattens
        s1 = srcrxn * NUM_SUFFIX
        s2 = srcrxn * REV_SUFFIX
        
        found = false
        
        # check for arm
        for arm_rxn in  [ARM_PREFFIX * srcrxn, ARM_PREFFIX * srcrxn * REV_SUFFIX]
            if arm_rxn in ec_rxns_set
                arm_rxndata = RxnData(ec_model, arm_rxn)
                push!(dict[srcrxn], arm_rxndata)
                delete!(ec_rxns_set, arm_rxn)
                found = true
            end
        end
        
        for ecrxn in ec_rxns_set
            if srcrxn == ecrxn || 
                startswith(ecrxn, s1) ||
                startswith(ecrxn, s2)
                
                # met data
                rxndata = RxnData(ec_model, ecrxn)
                push!(dict[srcrxn], rxndata)
                delete!(ec_rxns_set, ecrxn)
                
                # Extract draw and arm reaction
                for met in rxndata.mets
                    if startswith(met, PROT_PREFFIX) 
                        draw_rxn = DRAW_PREFFIX * met
                        draw_rxndata = RxnData(ec_model, draw_rxn)
                        push!(dict[srcrxn], draw_rxndata)
                        delete!(ec_rxns_set, draw_rxn)
                    end
                end
                
                # check for arm
                arm_rxn = ARM_PREFFIX * ecrxn
                if arm_rxn in ec_rxns_set
                    arm_rxndata = RxnData(ec_model, arm_rxn)
                    push!(dict[srcrxn], arm_rxndata)
                    delete!(ec_rxns_set, arm_rxn)
                end
                
                found = true
                                
            end
        end
        
        found && (count += 1)
        verbose && update!(prog, i) 

    end
    
    # Extras (Not related with any src reaction)
    dict[EXTRAS_KEY] = RxnData[]
    # prot_pool_exchange
    prot_pool_exchange = Human1.PROT_POOL_EXCHANGE
    if prot_pool_exchange in ec_rxns_set
        rxndata = RxnData(ec_model, prot_pool_exchange)
        push!(dict[EXTRAS_KEY], rxndata)
        delete!(ec_rxns_set, prot_pool_exchange)
    end
    
    verbose && finish!(prog) 
    return dict
end

# A reaction centric container
struct RxnData
    rxn::String
    subSys::Any
    mets::Vector{String}
    stoi::Vector{Float64}
    b::Vector{Float64}
    lb::Float64
    ub::Float64
    
    function RxnData(model::MetNet, ider)
        rxni = rxnindex(model, ider)
        rxn = model.rxns[rxni]
        subSys = model.subSystems[rxni]
        metis = rxn_mets(model, rxni)
        mets = model.mets[metis]
        stoi = [model.S[meti, rxni] for meti in metis]
        b = model.b[metis]
        lb, ub = model.lb[rxni], model.ub[rxni]
        return new(rxn, subSys, mets, stoi, b, lb, ub)
    end
end

get_ec_refdata(src_model, ec_model) = Dict(
        ECMAP_KEY => get_rxn_ecmap(src_model, ec_model; verbose = true),
        SRC_MODEL_KEY => src_model,
        EC_MODEL_KEY => ec_model,
        PROTLESS_KEY => ec_model.rxns[collect_protless(ec_model)],
        PROT_STOIS_KEY => prot_stois(ec_model)
    )