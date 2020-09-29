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