# This functions are intented to be use in conjuction for
# preparing a model for fba/EP. They are often not pure functions!

# redefining the Chemostat exchange criteria
function is_exchange_subsys(model, ider)
    idx = rxnindex(model, ider)
    subsys = model.subSystems[idx]
    return occursin(EXCH_SUBSYS_HINT, string(subsys))
end

"""
    delete_boundary_mets(base_model; verbose = true)

    I will delete all the Boundary (x) (see comps in the matmodel) metabilites, 
    leaving only the Extracellular (s) metabolites in the exchange reactions. 
    Why? they are not required
"""
function delete_boundary_mets(base_model; verbose = true)
    verbose && println("\nDeleting Boundary (x) metabolites ")
    verbose && println("Before: ", size(base_model))
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    base_model = del_met(base_model, to_del);
    verbose && println("After: ", size(base_model))
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    @assert isempty(to_del)
    return base_model
end

function prepare_extract_exchanges!(base_model::MetNet; verbose = true)
    # Exchanges
    exchs = []
    bkwd_exchs = []
    fwd_exchs = []

    subsys_exchs = filter((ider) -> is_exchange_subsys(base_model, ider), base_model.rxns)
    for exch in subsys_exchs
        exch_i = rxnindex(base_model, exch)

        # Because, this reactions are forward unbalanced (A <-> nothing)
        # positibe (+) bounds limit the outtake of the cell and
        # negative (-) bounds limit the intake.
        # Because in the Chemostat the intakes is 
        # controlled by the medium, we'll close all intakes to handle with them later
        # We'll open all outtakes
        react = rxn_reacts(base_model, exch_i)
        if isempty(react) 
            # backward defined (nothing -> A) 
            # A positive flux means intake (close)
            ub!(base_model, exch_i, 0.0)
            lb!(base_model, exch_i, -MAX_BOUND)
            push!(bkwd_exchs, exch)
        else 
            # forward defined (A -> nothing)
            # A positive flux means production (open)
            ub!(base_model, exch_i, MAX_BOUND)
            lb!(base_model, exch_i, 0.0)
            push!(fwd_exchs, exch)
        end

        push!(exchs, exch)
    end

    verbose && println("\nExchanges: ", exchs |> length)
    verbose && println("\tfwd_exchs (outputs): ", fwd_exchs |> length)
    verbose && println("\tbkwd_exchs (inputs): ", bkwd_exchs |> length)
    
    return exchs
end

function del_REV_rxns(model, rxns)
    println("\nDeleting REV exchs")
    println("Before: ", size(model))
    REV_rxns = []
    foreach(rxns) do rxn
        if endswith(rxn, REV_SUFFIX)
            push!(REV_rxns, rxn)
        else
            rxn *= REV_SUFFIX
            (rxn in model.rxns) && push!(REV_rxns, rxn)
        end
    end
    model = del_rxn(model, REV_rxns);
    println("After: ", size(model))
    prepare_extract_exchanges!(model)
    return model
end

function try_fba(args...; verbose = true)
    fbaout = nothing
    try
        fbaout = fba(args...)
        verbose && println("fba obj_val: ", fbaout.obj_val)
    catch err
        verbose && @warn("ERROR doing fba: ", err)
    end
    return fbaout
 end

 function apply_chstat_bound!(base_model, ξ, intake_info; verbose = true)
    verbose && println("\nApplying chemostat bound, xi: ", ξ)
    
    apply_bound!(base_model, ξ, intake_info; emptyfirst = true, ignore_miss = true);
    
    !verbose && return
    println("Medium exchanges: (-)intake/(+)output")
    for (rxn, lb) in base_model.intake_info
        eq = rxn_str(base_model, rxn);
        bs = bounds(base_model, rxn);
        println("\t", rxn, ": ", eq, " ", bs |> collect)
    end
end


