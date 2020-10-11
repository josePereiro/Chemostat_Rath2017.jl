# -*- coding: utf-8 -*-
using DrWatson
@quickactivate "Chemostat_Rath2017"

import DataFrames: DataFrame
import MAT
import CSV
import JSON

import Chemostat
const Ch = Chemostat
import Chemostat_Rath2017: HumanGEM, RathData, DATA_KEY
const HG = HumanGEM
const Rd = RathData

# This file is the primary input to the processing
if !isfile(HG.MODEL_RAW_MAT_FILE)
    error("$(relpath(HG.MODEL_RAW_MAT_FILE)) not found, you must run 'make all' fisrt (see README)!!!")
end

# ---
# ## Description
# ---

# This script will produce a file that contains the biomass reaction of the HumanGEM modified using the data extracted from Niklas (2013) (see comments)

# ---
# ### MAT model
# ---

mat_model = MAT.matread(HG.MODEL_RAW_MAT_FILE);
mat_model = mat_model[first(keys(mat_model))];
model = Ch.Utils.read_mat(HG.MODEL_RAW_MAT_FILE);
biomass_ider = "biomass_human"
biomass_idx = Ch.Utils.rxnindex(model, biomass_ider)
println("Loaded Model: ", size(model))

# ---
# ### Biomass equation
# ---
# I will modified the biomass equation (biomass_human) of Human1 model with data
# derived from Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1.
# I compute de relation between the total of each group reported in Niklas 2013 with the equivalent group found in the model biomass, and then rescaled each group to match the reported total. I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi

biomass = Dict()
for met_idx in Ch.Utils.rxn_mets(model, biomass_ider)
    met = model.mets[met_idx]
    biomass[met] = model.S[met_idx, biomass_idx]
end

# ### Carbohydrates

ch_ids = ["m03161c"]
exp_ch_tot = 438.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total ch: ", exp_ch_tot)
model_ch_tot = abs.(sum([biomass[met] for met in ch_ids])) # The model did not include directly any carbohydrate
println("model total ch: ", model_ch_tot)
ch_factor = exp_ch_tot/model_ch_tot
println("factor exp/model: ", ch_factor)

# ### RNA

rna_ids = ["m02847c"]
model_rna_tot = abs.(sum([biomass[met] for met in rna_ids]))
println("model total rna: ", model_rna_tot)
exp_rna_tot = 176.9 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total rna: ", exp_rna_tot)
rna_factor = exp_rna_tot/model_rna_tot
println("factor exp/model: ", rna_factor)

# ### DNA

dna_ids = ["m01721n"]
model_dna_tot = abs.(sum([biomass[met] for met in dna_ids]))
println("model total dna: ", model_dna_tot)
exp_dna_tot = 45.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total dna: ", exp_dna_tot)
dna_factor = exp_dna_tot/model_dna_tot
println("factor exp/model: ", dna_factor)

# ### Lipids

lip_ids = ["m10014c"]
model_lip_tot = abs.(sum([biomass[met] for met in lip_ids]))
println("model total lipids:", model_lip_tot)
# Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
exp_lip_tot = 202.9 * 1e-3
println("experimental total lipids:", exp_lip_tot)
lip_factor = exp_lip_tot/model_lip_tot
println("factor exp/model:", lip_factor)

# ### Aminoacids

aa_ids = ["m10013c"]
model_prot_tot = abs.(sum([biomass[met] for met in aa_ids]))
println("model total protein: ", model_prot_tot)
exp_prot_tot = 7462.6 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total protein: ", exp_prot_tot)
prot_factor = exp_prot_tot/model_prot_tot
println("factor exp/model: ", prot_factor)

# ### Rescaling

# +
# Carbohydrates 
for ch in ch_ids
    biomass[ch] = biomass[ch] * ch_factor
end

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

# Print
println("\nModified Biomass")
JSON.print(biomass, 4)
println()

# ### Saving

file = HG.NIKLAS_BIOMASS_FILE
tagsave(file, Dict(DATA_KEY => biomass))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


