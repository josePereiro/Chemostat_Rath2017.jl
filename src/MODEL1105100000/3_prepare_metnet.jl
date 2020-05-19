# TODO move this to prepare model
# ---
# ## Processing
# ---

# print()
# print('Processing')

# ### ATP maintinance demand


# # Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
# # pgDW/μm³ * 1e9 = pgDW/μL
# # pgDW/μL * 1e6 = pgDW/L
# # pgDW/L * 1e-12 = gDW/L
# # atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
# # # Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
# atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

# atpm_rxn = model.reactions.ATPmaint_LSQBKT_c_RSQBKT_
# print("ATP maintinance demand: ", atpm_rxn.id)
# print(atpm_rxn.reaction, atpm_rxn.bounds)
# print("lower bound set to: ", atpm_flux, "(mmol/gWD hr)")
# atpm_rxn.lower_bound = atpm_flux
# print(atpm_rxn.reaction, atpm_rxn.bounds)

# # ### Closing all intakes

# rxn = model.reactions.EX_2425dhvitd3_LPAREN_e_RPAREN_
# rxn.reaction



# for rxn in model.reactions:
#     if rxn.id.startswith("DM_"):
#         rxn.bounds = (0.0, 0.0)
#     if rxn.id.startswith("EX_"):
#         if rxn.reactants:
            
#         print(rxn.id)