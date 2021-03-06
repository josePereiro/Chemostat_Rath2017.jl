# +
import CSV
import DataFrames: DataFrame
import MAT

import Chemostat_Rath2017
H1 = Chemostat_Rath2017.Human1
# -

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# This create a map between the metabolites ids of the Human1 model 
# and the ids used in Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)

mets_map = Dict{AbstractString, AbstractString}()
mets_map["LAC"] = "m01716s"
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

# makinig it both ways
for (k, v) in mets_map
    mets_map[v] = k
end

# Saving
df = DataFrame([collect(keys(mets_map)), collect(values(mets_map))])
CSV.write(H1.METS_MAP_FILE, df)
println("created $(relpath(H1.METS_MAP_FILE))")



