# -*- coding: utf-8 -*-
import DataFrames: DataFrame
import MAT
import CSV
# import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
M = Chemostat_Rath2017.MODEL1105100000
Rd = Chemostat_Rath2017.RathData

# -
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# This file is the primary input to the processing
if !isfile(M.MODEL_RAW_MAT_FILE)
    error("$(M.MODEL_RAW_MAT_FILE) not found, you must run 'make all' fisrt (see README)!!!")
end

# ---
# ### MAT model
# ---

mat_model = MAT.matread(M.MODEL_RAW_MAT_FILE);
mat_model = mat_model[first(keys(mat_model))];
model = Ch.Utils.MetNet(mat_model);
println("Loaded Model: ", size(model))

biomass_ider = "BIOMASS"
biomass_idx = Ch.Utils.rxnindex(model, biomass_ider);



# ---
# ### Biomass equation
# ---
# I will modified the biomass equation of MODEL1105100000 model with data
# derived from Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1.
# I compute de relation between the total of each group reported in Niklas 2013 with the equivalent group found in the model biomass, and then rescaled each group to match the reported total. I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi

biomass = Dict()
for met_idx in Ch.Utils.rxn_mets(model, biomass_ider)
    met = model.mets[met_idx]
    biomass[met] = model.S[met_idx, biomass_idx]
end

# ### Carbohydrates

ch_id = "g6p_c"
exp_ch_tot = 438.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total ch: ", exp_ch_tot)
model_ch_tot = 0.0 # The model did not include directly any carbohydrate
println("model total ch: ", model_ch_tot)

# ### RNA

rna_ids = ["amp_c", "cmp_c", "gmp_c", "ump_c"]
model_rna_tot = abs.(sum([biomass[met] for met in rna_ids]))
println("model total rna: ", model_rna_tot)
exp_rna_tot = 176.9 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total rna: ", exp_rna_tot)
rna_factor = exp_rna_tot/model_rna_tot
println("factor exp/model: ", rna_factor)

# ### DNA

dna_ids = ["damp_c", "dcmp_c", "dgmp_c", "dtmp_c"]
model_dna_tot = abs.(sum([biomass[met] for met in dna_ids]))
println("model total dna: ", model_dna_tot)
exp_dna_tot = 45.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total dna: ", exp_dna_tot)
dna_factor = exp_dna_tot/model_dna_tot
println("factor exp/model: ", dna_factor)

# ### Lipids

lip_ids = ["chsterol_c", "clpn_DASH_hs_c", "dag_DASH_hs_c",
           "lpchol_DASH_hs_c", "mag_DASH_hs_c", "pa_DASH_hs_c",
           "pail_DASH_hs_c", "pe_DASH_hs_c", "ps_DASH_hs_c",
           "sphmyln_DASH_hs_c", "tag_DASH_hs_c", "xolest_DASH_hs_c"]
model_lip_tot = abs.(sum([biomass[met] for met in lip_ids]))
println("model total lipids:", model_lip_tot)
# Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
exp_lip_tot = 202.9 * 1e-3
println("experimental total lipids:", exp_lip_tot)
lip_factor = exp_lip_tot/model_lip_tot
println("factor exp/model:", lip_factor)

# ### Aminoacids

aa_ids = ["ala_DASH_DASH_DASH_L_c", "asn_DASH_DASH_DASH_L_c", "asp_DASH_DASH_DASH_L_c",
          "gln_DASH_DASH_DASH_L_c", "glu_DASH_DASH_DASH_L_c", "pro_DASH_DASH_DASH_L_c",
          "ser_DASH_DASH_DASH_L_c", "gly_c"]
model_prot_tot = abs.(sum([biomass[met] for met in aa_ids]))
println("model total protein: ", model_prot_tot)
exp_prot_tot = 7462.6 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total protein: ", exp_prot_tot)
prot_factor = exp_prot_tot/model_prot_tot
println("factor exp/model: ", prot_factor)

# ### Rescaling

# +
# Carbohydrates # I put all the reported quatity in a single metabolite
biomass[ch_id] = -exp_ch_tot

# Aminoacids
for aa in aa_ids
    biomass[aa] = biomass[aa] * prot_factor
end

# lipids
for lip in lip_ids
    biomass[lip] = biomass[lip] * lip_factor
end

# DNA
for dna in dna_ids
    biomass[dna] = biomass[dna] * dna_factor
end

# RNA
for rna in rna_ids
    biomass[rna] = biomass[rna] * rna_factor
end
# -

# ### Saving

df = DataFrame([collect(keys(biomass)), collect(values(biomass))])
CSV.write(M.NIKLAS_BIOMASS_FILE, df)
println("created $(relpath(M.NIKLAS_BIOMASS_FILE))")


