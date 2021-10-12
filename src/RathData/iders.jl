# TODO: Make uppercase
# ids related with experiments performed by Rath2017
const exps = ["A", "B", "C", "D", "E", "F01","F02","F03"] 
# experiments that reached steady states, ordered by ξ value
const ststs = ["E", "D", "F01", "A", "B", "C"]
# measured mets
# Mets which concentration and exchange was measured
const msd_mets = ["GLC", "LAC", "GLN", "NH4", "GAL","PYR","GLU","ALA","ASP"]
# all_mets
# Counting the feed medium mets and the measured mets
const all_mets = ["LAC", "NH4", "HCO3", "ASN", "MET",
    "GLC", "GLY", "TYR", "LEU", "PHE", "GLU", "ASP",
    "GLN", "CYS", "ILE", "PYR", "H", "LYS", "TRP", "SER",
    "HIS", "VAL", "ALA", "ARG", "PRO", "GAL", "THR"]
# exchange rates
const rate_q = ["q$met" for met in msd_mets]
# feed concentrations
const feed_c = ["cGLC", "cGAL", "cGLN"]
# medium concentration
const med_s = ["s$met" for met in msd_mets];
const growth_ider = "μ"
const iders_to_plot = [msd_mets; growth_ider]