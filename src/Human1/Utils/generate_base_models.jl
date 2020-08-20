# This functions are intented to be use in conjuction for
# preparing a model for fba/EP. They are often not pure functions!

# redefining the Chemostat exchange criteria
function is_exchange(model, ider, 
    exch_subsys_hint = "Exchange/demand")
    idx = rxnindex(model, ider)
    subsys = model.subSystems[idx]
    return occursin(exch_subsys_hint, string(subsys))
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

function prepare_extract_exchanges!(base_model::MetNet, max_bound = MAX_BOUND; verbose = true)
    # Exchanges
    exchs = []
    bkwd_exchs = []
    fwd_exchs = []
    for exch in filter((ider) -> is_exchange(base_model, ider), base_model.rxns)
        exch_i = rxnindex(base_model, exch)

        # First close it, later if it is what I want, open the outtake
        lb!(base_model, exch_i, 0.0) 
        ub!(base_model, exch_i, 0.0)

        mets = rxn_mets(base_model, exch_i)

        length(mets) != 1 && continue # I want only monomoleculars
        !endswith(base_model.mets[first(mets)], "s") && continue # I what only the exchanges 's'

        # Because, this reactions are forward unbalanced (A <-> nothing)
        # positibe (+) bounds limit the outtake of the cell and
        # negative (-) bounds limit the intake.
        # Because in the Chemostat the intakes is 
        # controlled by the medium, we'll close all intakes to handle with them later
        # We'll open all outtakes
        react = rxn_reacts(base_model, exch_i)
        if isempty(react) 
            # backward defined (nothing <- A) 
            # A positive flux means intake (close)
            push!(bkwd_exchs, exch_i)
            ub!(base_model, exch_i, 0.0)
        else # forward defined (A -> nothing)
            # A positive flux means production (open)
            push!(fwd_exchs, exch_i)
            ub!(base_model, exch_i, max_bound)
        end

        push!(exchs, exch_i)
    end
    verbose && println("\nExchanges: ", exchs |> length)
    verbose && println("\tfwd_exchs (outputs): ", fwd_exchs |> length)
    verbose && println("\tbkwd_exchs (inputs): ", bkwd_exchs |> length)
    
    return (exchs, bkwd_exchs, fwd_exchs)
end

function del_bkwd_exchs(base_model, bkwd_exchs, max_bound = MAX_BOUND)
    println("\nDeleting bkwd_exchs")
    println("Before: ", size(base_model))
    base_model = del_rxn(base_model, bkwd_exchs);
    println("After: ", size(base_model))
    exchs, bkwd_exchs, fwd_exchs = prepare_extract_exchanges!(base_model, max_bound)
    return base_model, exchs, bkwd_exchs, fwd_exchs
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
    verbose && println("\nApplaying chemostat bound, xi: ", ξ)
    
    apply_bound!(base_model, ξ, intake_info; emptyfirst = true, ignore_miss = true);
    
    !verbose && return
    println("Medium exchanges: (-)intake/(+)output")
    for (rxn, lb) in base_model.intake_info
        eq = rxn_str(base_model, rxn);
        bs = bounds(base_model, rxn);
        println("\t", rxn, ": ", eq, " ", bs |> collect)
    end
end
