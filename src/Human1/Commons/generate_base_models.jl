# This functions are intented to be use in conjuction for
# preparing a model for fba/EP. They are often not pure functions!

# redefining the Chemostat exchange criteria
function is_exchange_subsys(model, ider)
    idx = ChU.rxnindex(model, ider)
    subsys = model.subSystems[idx]
    return occursin(EXCH_SUBSYS_HINT, string(subsys))
end

subsys_exchs(model) = filter((ider) -> is_exchange_subsys(model, ider), eachindex(model.rxns))

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
    base_model = ChU.del_met!(base_model, to_del);
    base_model = ChU.compacted_model(base_model)
    verbose && println("After: ", size(base_model))
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    @assert isempty(to_del)
    return base_model
end

function prepare_extract_exchanges!(base_model::ChU.MetNet; verbose = true)
    # Exchanges
    exchs = []
    bkwd_exchs = []
    fwd_exchs = []

    ss_exchs = base_model.rxns[subsys_exchs(base_model)]
    for exch in ss_exchs
        exch_i = ChU.rxnindex(base_model, exch)

        # Because, this reactions are forward unbalanced (A <-> nothing)
        # positibe (+) bounds limit the outtake of the cell and
        # negative (-) bounds limit the intake.
        # Because in the Chemostat the intakes is 
        # controlled by the medium, we'll close all intakes to handle with them later
        # We'll open all outtakes
        react = ChU.rxn_reacts(base_model, exch_i)
        if isempty(react) 
            # backward defined (nothing -> A) 
            # A positive flux means intake (close)
            ChU.ub!(base_model, exch_i, 0.0)
            ChU.lb!(base_model, exch_i, -MAX_BOUND)
            push!(bkwd_exchs, exch)
        else 
            # forward defined (A -> nothing)
            # A positive flux means production (open)
            ChU.ub!(base_model, exch_i, MAX_BOUND)
            ChU.lb!(base_model, exch_i, 0.0)
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
    ChU.del_rxn!(model, REV_rxns)
    model = ChU.compacted_model(model)
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

function apply_biomass!(base_model, obj_ider, niklas_biomass)
    # I will modified the biomass equation of MODEL1105100000 model with data
    # derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
    # I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
    println("\nApplying Niklas Biomass")
    # JSON.print(HG.niklas_biomass, 4)
    # println()

    biomass_idx = rxnindex(base_model, obj_ider)
    base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
    for (met, y) in niklas_biomass
        S!(base_model, met, biomass_idx, y)
    end
end

function set_atp_demand(base_model)
    base_model = deepcopy(base_model);
    # From HumanGEM ATPM
    # SUMMARY (color code: warning, info, error)
    #  HMR_3964 ()
    #  lb: 0.0, ub: 1000.0
    #  (-1.0) m01371c + (-1.0) m02040c ==> (1.0) m01285c + (1.0) m02039c + (1.0) m02751c
    atpm_ider = Human1.ATPM_IDER
    if !(atpm_ider in base_model.rxns)
        base_model = add_rxn(base_model, atpm_ider; 
                mets = Dict("m01371c" => -1.0, "m02040c" => -1.0, # ->
                    "m01285c" => 1.0, "m02039c" => 1.0, "m02751c" => 1.0), 
                lb = -MAX_BOUND,
                ub = MAX_BOUND);
    end
    # This reaction is foward defined with respect to atp
    # atp + more_reacts <-> adp + more_products
    # so we limit the lower bounds as the minimum atp demand 
    # that the cell metabolism must fullfill
    lb!(base_model, atpm_ider, ATPM_FLUX)

    println("\nATP demand")
    summary(base_model, atpm_ider)

    return base_model
end

function open_rxns!(model, iders)
    for rxn in iders
        bounds!(model, rxn, -MAX_BOUND, MAX_BOUND);
    end
end

# function run_fba_test(base_model, obj_ider)
#     fbaout = Ch.LP.fba(base_model, obj_ider);
#     println("\nFBAout summary")
#     Ch.Utils.summary(base_model, fbaout)

#     println("\nComparing with experiments")
#     model = deepcopy(base_model)
#     for stst in Rd.ststs
#         println("\nStst: $stst")
#         ξ = Rd.val(:ξ, stst)
#         println("exp xi: $ξ")
#         exp_μ = Rd.val(:μ, stst)
#         println("exp growth: $exp_μ")    
#         Ch.SteadyState.apply_bound!(model, ξ, Dict());
#         fbaout_μ = Ch.LP.fba(model, obj_ider).obj_val;
#         fbaout_μ < exp_μ && @warn("fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)")
#         println("fba growth: $fbaout_μ")
#     end
# end