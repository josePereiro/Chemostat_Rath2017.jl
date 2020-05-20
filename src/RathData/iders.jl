# ids related with experiments performed by Rath2017
exps = ["A", "B", "C", "D", "E", "F01","F02","F03"] 
# experiments that reached steady states, ordered by ξ value
ststs = ["E", "D", "F01", "A", "B", "C"]
# measured mets
# Mets which concentration and exchange was measured
msd_mets = ["GLC", "LAC", "GLN", "NH4", "GAL","PYR","GLU","ALA","ASP"]
# all_mets
# Counting the feed medium mets and the measured mets
all_mets = ["LAC", "NH4", "HCO3", "ASN", "MET",
    "GLC", "GLY", "TYR", "LEU", "PHE", "GLU", "ASP",
    "GLN", "CYS", "ILE", "PYR", "H", "LYS", "TRP", "SER",
    "HIS", "VAL", "ALA", "ARG", "PRO", "GAL", "THR"]
# exchange rates
rate_q = ["q$met" for met in msd_mets]
# feed concentrations
feed_c = ["cGLC", "cGAL", "cGLN"]
# medium concentration
med_s = ["s$met" for met in msd_mets];