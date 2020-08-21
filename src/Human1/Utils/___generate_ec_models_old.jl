"""
    Store all the analysis relevant data
"""
struct EcMGDC # EcModel Generation Data Container
    orig_model::MetNet
    ec_template::MetNet
    orig_rxns::Vector{Tuple{Int64,String}}
    orig_mets::Vector{Tuple{Int64,String}}
    ec_rxns::Vector{Tuple{Int64,String}}
    ec_mets::Vector{Tuple{Int64,String}}
    orig_rxns_left::Dict{Int64, String}
    orig_mets_left::Dict{Int64, String}
    ec_rxns_left::Dict{Int64, String}
    ec_mets_left::Dict{Int64, String}
    
    function EcMGDC(orig_model::MetNet, ec_template::MetNet)
        # Get cost related reactions
        # Here I'll keep all the data to be included in the ec new model
        orig_rxns = Tuple{Int64,String}[]
        orig_mets = Tuple{Int64,String}[]
        ec_rxns = Tuple{Int64,String}[]
        ec_mets = Tuple{Int64,String}[]

        # To track not included
        orig_rxns_left = Dict(i => rxn for (i, rxn) in orig_model.rxns |> enumerate) 
        orig_mets_left = Dict(i => met for (i, met) in orig_model.mets |> enumerate)
        ec_rxns_left = Dict(i => rxn for (i, rxn) in ec_template.rxns |> enumerate) 
        ec_mets_left = Dict(i => met for (i, met) in ec_template.mets |> enumerate)
        
        new(orig_model, ec_template,
            orig_rxns, orig_mets, ec_rxns, ec_mets, 
            orig_rxns_left, orig_mets_left, ec_rxns_left, ec_mets_left)
    end
end

_gd_str(gd::EcMGDC) = string("in/left ", 
                    "rxns[", 
                        "ec:",      length(gd.ec_rxns), 
                        "/",  length(gd.ec_rxns_left) , 
                        " orig:" , length(gd.orig_rxns),
                        "/" ,length(gd.orig_rxns_left),
                    "] mets[",
                        "ec:", length(gd.ec_mets), 
                        "/", length(gd.ec_mets_left), 
                        " orig:" , length(gd.orig_mets),
                        "/" , length(gd.orig_mets_left),
                    "]"
)
Base.show(io::IO, gd::EcMGDC) = print(io, _gd_str(gd));


function collect_rxns_data!(gd::EcMGDC; verbose = true, 
        # This controls if an original rxn must be 
        # search in the ec reference model
        on_skip_ec = function(gd, oi) 
            subsys = gd.orig_model.subSystems[oi];
            return occursin(EXCH_SUBSYS_HINT, string(subsys)) 
        end
    )
    verbose && println("Collecting rxns data")

    olen = length(gd.orig_model.rxns)
    for (oi, orxn) in gd.orig_model.rxns |> enumerate

        s1 = orxn * "No"
        s2 = orxn * "_REV"
        s3 = "arm_" * orxn

        found_in_ec = false       
        ec_dat = on_skip_ec(gd, oi) ? 
                empty(gd.ec_rxns_left) : 
                gd.ec_rxns_left
        for (eci, ecrxn) in ec_dat
            if orxn == ecrxn || 
                startswith(ecrxn, s1) ||
                startswith(ecrxn, s2) ||
                startswith(ecrxn, s3)

                found_in_ec = true

                push!(gd.ec_rxns, (eci, ecrxn))
                delete!(gd.ec_rxns_left, eci)

                # Just printing progress
                verbose && mod(oi, UP_FREC) == 0 && 
                    (Core.print("[", oi, " / ", olen, "] ", _gd_str(gd),
                        " ex: ", orxn, " => ", ecrxn, " "^20, "\r"); flush(stdout))

            end
        end

        if !found_in_ec
            push!(gd.orig_rxns, (oi, orxn))
            delete!(gd.orig_rxns_left, oi)
        end
    end

    unique!(gd.orig_rxns)
    unique!(gd.ec_rxns)

    verbose && println(" "^50, "\r\tDone: ", _gd_str(gd), " "^50)
end

function collect_mets_data!(gd::EcMGDC; verbose = true)
    verbose && println("Collecting mets data")
    
    for (rxnsdat, metsdat, leftmet, model) in [
                                    (gd.orig_rxns, gd.orig_mets, gd.orig_mets_left, gd.orig_model),
                                    (gd.ec_rxns  , gd.ec_mets  , gd.ec_mets_left  , gd.ec_template  )
                                ]

        len = length(rxnsdat)
        for (i, (rxni, rxn)) in rxnsdat |> enumerate

            metidxs = rxn_mets(model, rxni)
            metiders = model.mets[metidxs]
            foreach(zip(metidxs, metiders)) do (metidx, metider)
                push!(metsdat, (metidx, metider))
                delete!(leftmet, metidx)
            end


            if verbose && mod(i, UP_FREC) == 0
                unique!(metsdat)
                 # Just printing progress
                Core.print("[", i, " / ", len, "] ", _gd_str(gd), "\r")
                flush(stdout)
            end
        end
        unique!(metsdat)
    end

    verbose && println(" "^50, "\r\tDone: ", _gd_str(gd), " "^50)
end

function collect_draw_rxns!(gd::EcMGDC; verbose = true)
    # get involved prots pseudo-metabolites
    verbose && println("Collecting draw rxns prot_pool -> prot_X")
    prot_mets = filter(gd.ec_mets) do (meti, met)
        startswith(met, "prot_") && met != PROT_POOL
    end
    
    # get draw rxns prot_pool -> prot_X
    foreach(prot_mets) do (prot_meti, prot_met)
        draw_rxn = "draw_" * prot_met
        draw_rxni = rxnindex(gd.ec_template, draw_rxn)
        delete!(gd.ec_rxns_left, draw_rxni)
        push!(gd.ec_rxns, (draw_rxni, draw_rxn))
    end
    verbose && println("\tDone: ", _gd_str(gd), " "^50)
end

function add_prot_pool_exchange!(gd::EcMGDC; verbose = true)
    # Add prot_pool_exchange
    verbose && println("Add prot_pool_exchange")
    
    rxni = rxnindex(gd.ec_template, PROT_POOL_EXCHANGE)
    delete!(gd.ec_rxns_left, rxni)
    push!(gd.ec_rxns, (rxni, PROT_POOL_EXCHANGE))

    meti = metindex(gd.ec_template, PROT_POOL)
    delete!(gd.ec_mets_left, meti)
    push!(gd.ec_mets, (meti, PROT_POOL))
    verbose && println("\tDone: ", _gd_str(gd), " "^50)
end

function make_all_unique!(gd::EcMGDC; verbose = true)
    verbose && println("Make all unique (just for be sure)")
    
    verbose && println("\tBefore: ", _gd_str(gd), " "^50)
    for (ec_dat, orig_dat) in [(gd.ec_rxns, gd.orig_rxns), 
                               (gd.ec_mets, gd.orig_mets)]
        
        unique!(ec_dat)
        unique!(orig_dat)
        ec_iders = [ider for (i, ider) in ec_dat]
        orig_dat_ = copy(orig_dat)
        empty!(orig_dat)
        for (oi, oider) in orig_dat_
            oider in ec_iders && continue
            push!(orig_dat, (oi, oider))
        end
    end
    verbose && println("\tAfter: ", _gd_str(gd), " "^50)
end

# TODO: build a more complite model
function build_new_model(gd::EcMGDC)
    
    M, N = (length(gd.orig_mets) + length(gd.ec_mets), length(gd.orig_rxns) + length(gd.ec_rxns));

    S_ = zeros(M, N);
    b_ = zeros(M);
    lb_ = zeros(N);
    ub_ = zeros(N);
    rxns_ = Vector{String}(undef, N);
    mets_ = Vector{String}(undef, M);
    subSystems_ = Vector(undef, N);

    # This map between ider and new model idx
    rxns_newi_map = Dict(rxn => newi for (newi, (modi, rxn)) in 
                                    [gd.orig_rxns; gd.ec_rxns] |> enumerate)
    mets_newi_map = Dict(met => newi for (newi, (modi, met)) in 
                                    [gd.orig_mets; gd.ec_mets] |> enumerate);
    
    # Adding mets
    for (metsdat, model) in [(gd.orig_mets, gd.orig_model), 
                             (gd.ec_mets, gd.ec_template)]
        for (modmi, met) in metsdat
            newmi = mets_newi_map[met]

            # include met data
            mets_[newmi] = met
            b_[newmi] = b(model, modmi)
        end
    end
    
    # Adding rxns
    for (rxnsdat, model) in [(gd.orig_rxns, gd.orig_model), 
                            (gd.ec_rxns, gd.ec_template)]

        for (modri, rxn) in rxnsdat
            newri = rxns_newi_map[rxn]

            # include rxn data (including the stoichiometry)
            rxns_[newri] = rxn
            subSystems_[newri] = model.subSystems[modri]
            lb_[newri], ub_[newri] = bounds(model, modri)
            metis = rxn_mets(model, modri)
            for modmi in metis
                met = model.mets[modmi]
                newmi = mets_newi_map[met]
                S_[newmi, newri] = S(model, modmi, modri)
            end
        end
    end
    
    new_ec_model = MetNet(S_, b_, lb_, ub_, rxns_, mets_; subSystems = subSystems_);
end

function collect_protless(model::MetNet, idxs = eachindex(model.rxns))
    filter(idxs) do rxni
        metis = rxn_mets(model, rxni)
        metids = model.mets[metis]
        !any(startswith(met, "prot_") for met in metids)
    end
end

function prot_stois(ref_model)
    prot_kin_stois = []
    prot_draw_stois = []

    prot_pool_idx = metindex(ref_model, PROT_POOL)
    foreach(eachindex(ref_model.mets)) do meti
        met = ref_model.mets[meti]
        !startswith(met, "prot_") && return
        meti == prot_pool_idx && return

        rxnis = met_rxns(ref_model, meti)
        for rxni in rxnis
            rxn = ref_model.rxns[rxni]
            if startswith(rxn, "draw_")
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

