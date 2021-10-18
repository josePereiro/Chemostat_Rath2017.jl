
## ---------------------------------------------------------------------
# Base intake info
# The Base model will have a medium (expressed as open intake fluxes) that resamble the cultivation
# at xi = 1 using the set up in Rath 2017 exp A. So the lb of the intakes will be directly the (negative)
# concentration in 42_MAX_UB standard medium (see Cossio's paper). Also, a few intakes, not justified in
# the standard medium will be add based in the intakes of the original model FBA analysis.
function load_base_intake_info(; inf_medium = true)

    base_intake_info = Dict()
    mets_map = load_mets_map()
    met_readable_ids = load_met_readable_ids()
    exch_met_map = load_exch_met_map()

    # From 42_MAX_UB standard medium
    for rath_met in RathData.ALL_METS
        
        # base_model id
        model_met = mets_map[rath_met]   
        exch_rxn = exch_met_map[model_met]
        
        # 42_MAX_UB standard medium
        # we take the openest version of the intakes for building the
        # base model
        conc = maximum(RathData.cval(rath_met, RathData.STSTS, 0.0)) 
        lb = -ABS_MAX_BOUND # intake bound
        base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
    end

    # From ham's medium see Hams_medium.jl
    for (Ham_id, conc) in load_Ham_medium(;inf_medium = inf_medium)
        model_met = met_readable_ids[Ham_id]
        exch_rxn = exch_met_map[model_met]
        haskey(base_intake_info, exch_rxn) && continue # Not overwrite 42_MAX_UB standard medium
        
        lb = -ABS_MAX_BOUND # intake bound
        base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
    end

    return base_intake_info

end

"""
    returns a copy of the base_intake_info but with the 
    feed medium concentration of a given Rath steady state
"""
function stst_base_intake_info(stst) 
    intake_info = load_base_intake_info()
    mets_map = load_mets_map()
    exch_met_map = load_exch_met_map()
    
    # The feed medium of each steady state only vary
    # in the composition of this mets (see Rath2017)
    for rath_met in ["GLC", "GLN", "GAL"]
        model_met = mets_map[rath_met]
        model_exch = exch_met_map[model_met]
        intake_info[model_exch]["c"] = RathData.cval(rath_met, stst)
    end
    return intake_info
end