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
        intake_info[model_exch]["c"] = RathData.val("c$rath_met", stst)
    end
    @assert intake_info["EX_glc_LPAREN_e_RPAREN"]["c"] == RathData.val("cGLC", stst)
    return intake_info
end

# TODO: move to scripts
# ### add demand
# function add_a1at_synthesis!(metnet, stst)
#     for (aa, rel_ab) in RathData.a1at_aa_rel_ab
#         model_id = mets_map[aa][1:end-1] * "c" # Change to internal metabolite
#         # We negate the flux becouse it must be a possitive demand
#         e = rel_ab * -RathData.val("qA1AT", stst)
#         Ch.Utils.b!(metnet, model_id, e)
#     end
# end