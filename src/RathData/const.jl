# TODO: Make uppercase
# ids related with experiments performed by Rath2017
const EXPS = ["A", "B", "C", "D", "E", "F01","F02","F03"] 
# experiments that reached steady states, ordered by ξ value
const STSTS = ["E", "D", "F01", "A", "B", "C"]
# measured mets
# Mets which concentration and exchange was measured
const MSD_METS = ["GLC", "LAC", "GLN", "NH4", "GAL","PYR","GLU","ALA","ASP"]
# ALL_METS
# Counting the feed medium mets and the measured mets
const ALL_METS = ["LAC", "NH4", "HCO3", "ASN", "MET",
    "GLC", "GLY", "TYR", "LEU", "PHE", "GLU", "ASP",
    "GLN", "CYS", "ILE", "PYR", "H", "LYS", "TRP", "SER",
    "HIS", "VAL", "ALA", "ARG", "PRO", "GAL", "THR"]
# exchange rates
const RATE_Q = ["q$met" for met in MSD_METS]
# feed concentrations
const FEED_C = ["cGLC", "cGAL", "cGLN"]
# medium concentration
const MED_S = ["s$met" for met in MSD_METS];
const GROWTH_IDER = "μ"