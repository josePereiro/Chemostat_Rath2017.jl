# +
using DrWatson
@quickactivate "Chemostat_Rath2017"

import JSON # for pretty printing

import Chemostat_Rath2017: DATA_KEY, HumanGEM
const HG = HumanGEM
# -

# This create a map between the metabolites ids of the Human1 model 
# and the ids used in Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)

mets_map = Dict{AbstractString, AbstractString}()
# Mammalian cells only contain L-LDH so that in humans the lactate produced is almost exclusively L-lactate.
# https://acutecaretesting.org/en/articles/l-lactate-and-d-lactate-clinical-significance-of-the-difference
mets_map["LAC"] = "m02403s" # "L-lactate[s]" # "m01716s" # "D-lactate[s]"
mets_map["NH4"] = "m02579s"
mets_map["HCO3"] = "m02046s"
mets_map["ASN"] = "m01369s"
mets_map["MET"] = "m02471s"
mets_map["GLC"] = "m01965s"
mets_map["GLY"] = "m01986s"
mets_map["TYR"] = "m03101s"
mets_map["LEU"] = "m02360s"
mets_map["PHE"] = "m02724s"
mets_map["GLU"] = "m01974s"
mets_map["ASP"] = "m01369s"
mets_map["GLN"] = "m01975s"
mets_map["CYS"] = "m01628s"
mets_map["ILE"] = "m02184s"
mets_map["PYR"] = "m02819s"
mets_map["H"] = "m02039s"
mets_map["LYS"] = "m02426s"
mets_map["TRP"] = "m03089s"
mets_map["SER"] = "m02896s"
mets_map["HIS"] = "m02125s"
mets_map["VAL"] = "m03135s"
mets_map["ALA"] = "m01307s"
mets_map["ARG"] = "m01365s"
mets_map["PRO"] = "m02770s"
mets_map["GAL"] = "m01910s"
mets_map["THR"] = "m02993s"

# print
println("\nRath2017 - HumanGEM met map")
JSON.print(mets_map, 4)
println()

# makinig it both ways
for (k, v) in mets_map
    mets_map[v] = k
end

# Saving
file = HG.METS_MAP_FILE
tagsave(file, Dict(DATA_KEY => mets_map))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")



