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