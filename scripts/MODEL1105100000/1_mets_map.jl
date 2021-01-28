import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

import UtilsJL
const UJL = UtilsJL

import Chemostat_Rath2017
const ChR = Chemostat_Rath2017
const M = ChR.MODEL1105100000

## ----------------------------------------------------------------------------
# This create a map between the metabolites ids of the MODEL1105100000_url model 
# and the ids used in Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)

mets_map = Dict{AbstractString, AbstractString}()
mets_map["LAC"] = "lac_DASH_DASH_DASH_L_e"
mets_map["NH4"] = "nh4_e"
mets_map["HCO3"] = "hco3_e"
mets_map["ASN"] = "asn_DASH_DASH_DASH_L_e"
mets_map["MET"] = "met_DASH_DASH_DASH_L_e"
mets_map["GLC"] = "glc_DASH_DASH_DASH_D_e"
mets_map["GLY"] = "gly_e"
mets_map["TYR"] = "tyr_DASH_DASH_DASH_L_e"
mets_map["LEU"] = "leu_DASH_DASH_DASH_L_e"
mets_map["PHE"] = "phe_DASH_DASH_DASH_L_e"
mets_map["GLU"] = "glu_DASH_DASH_DASH_L_e"
mets_map["ASP"] = "asp_DASH_DASH_DASH_L_e"
mets_map["GLN"] = "gln_DASH_DASH_DASH_L_e"
mets_map["CYS"] = "cys_DASH_DASH_DASH_L_e"
mets_map["ILE"] = "ile_DASH_DASH_DASH_L_e"
mets_map["PYR"] = "pyr_e"
mets_map["H"] = "h_e"
mets_map["LYS"] = "lys_DASH_DASH_DASH_L_e"
mets_map["TRP"] = "trp_DASH_DASH_DASH_L_e"
mets_map["SER"] = "ser_DASH_DASH_DASH_L_e"
mets_map["HIS"] = "his_DASH_DASH_DASH_L_e"
mets_map["VAL"] = "val_DASH_DASH_DASH_L_e"
mets_map["ALA"] = "ala_DASH_DASH_DASH_L_e"
mets_map["ARG"] = "arg_DASH_DASH_DASH_L_e"
mets_map["PRO"] = "pro_DASH_DASH_DASH_L_e"
mets_map["GAL"] = "gal_e"
mets_map["THR"] = "thr_DASH_DASH_DASH_L_e"

# makinig it both ways
for (k, v) in mets_map
    mets_map[v] = k
end

## ----------------------------------------------------------------------------
# Saving
UJL.save_data(M.METS_MAP_FILE, mets_map)

