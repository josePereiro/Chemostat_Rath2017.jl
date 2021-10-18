using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin
    import DataFrames: DataFrame
    import MAT
    import CSV

    import Chemostat
    import Chemostat.MetNets
    const Ch = Chemostat
    import Chemostat_Rath2017: Human1, RathData
    const Rd = RathData
    const H1 = Human1
    const HG = H1.HumanGEM
end

## ------------------------------------------------------------------
# load model
model = HG.load_humangem_raw_model()
biomass_idx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER)

## ------------------------------------------------------------------
# Biomass equation

# I will modified the biomass equation (biomass_human) of Human1 model with data
# derived from Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1.
# I compute de relation between the total of each group reported in Niklas 2013 with the equivalent
# group found in the model biomass, and then rescaled each group to match the reported total.
# I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi

biomass = Dict()
for met_idx in MetNets.rxn_mets(model, biomass_idx)
    met = model.mets[met_idx]
    biomass[met] = model.S[met_idx, biomass_idx]
end

## ------------------------------------------------------------------
# Carbohydrates

ch_ids = ["m03161c"]
exp_ch_tot = 438.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total ch: ", exp_ch_tot)
model_ch_tot = abs.(sum([biomass[met] for met in ch_ids])) # The model did not include directly any carbohydrate
println("model total ch: ", model_ch_tot)
ch_factor = exp_ch_tot/model_ch_tot
println("factor exp/model: ", ch_factor)

## ------------------------------------------------------------------
# RNA

rna_ids = ["m02847c"]
model_rna_tot = abs.(sum([biomass[met] for met in rna_ids]))
println("model total rna: ", model_rna_tot)
exp_rna_tot = 176.9 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total rna: ", exp_rna_tot)
rna_factor = exp_rna_tot/model_rna_tot
println("factor exp/model: ", rna_factor)

## ------------------------------------------------------------------
# DNA

dna_ids = ["m01721n"]
model_dna_tot = abs.(sum([biomass[met] for met in dna_ids]))
println("model total dna: ", model_dna_tot)
exp_dna_tot = 45.3 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total dna: ", exp_dna_tot)
dna_factor = exp_dna_tot/model_dna_tot
println("factor exp/model: ", dna_factor)

## ------------------------------------------------------------------
# Lipids

lip_ids = ["m10014c"]
model_lip_tot = abs.(sum([biomass[met] for met in lip_ids]))
println("model total lipids:", model_lip_tot)
# Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
exp_lip_tot = 202.9 * 1e-3
println("experimental total lipids:", exp_lip_tot)
lip_factor = exp_lip_tot/model_lip_tot
println("factor exp/model:", lip_factor)

# Aminoacids

aa_ids = ["m10013c"]
model_prot_tot = abs.(sum([biomass[met] for met in aa_ids]))
println("model total protein: ", model_prot_tot)
exp_prot_tot = 7462.6 * 1e-3 # Niklas (2013): 103–114. https://doi.org/10.1016/j.ymben.2013.01.002. Table1
println("experimental total protein: ", exp_prot_tot)
prot_factor = exp_prot_tot/model_prot_tot
println("factor exp/model: ", prot_factor)

## ------------------------------------------------------------------
# Rescaling

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

## ------------------------------------------------------------------
# save
sdat(HG, biomass,
    "niklas_biomass", ".jls";
    verbose = true
)