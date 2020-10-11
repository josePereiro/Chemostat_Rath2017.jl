## ---------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## ---------------------------------------------------------------
import JSON # For pretty printing

import Chemostat_Rath2017: Chemostat, Human1
const ChU = Chemostat.Utils
const HG = Human1.HumanGEM


## ---------------------------------------------------------------
# Ham's F-12 
# but mediam componets from Human1 Zenodo https://zenodo.org/record/3583004#.Xu8KdmpKiCh
# task (Growth in Ham's medium)
# Final concentrations in mM
Inf_ = 9999 # taked as infinite concentration
def_ = 1e-1 #

# Nutrient Mixture F-12 Ham N6760
# Source concentration where are g/L, from https://www.sigmaaldrich.com/life-science/cell-culture/learning-center/media-formulations/f-12-ham.html

Ham_medium = Dict()

## ---------------------------------------------------------------
# Inorganic Salts

# Calcium Chloride (CaCl2)   0.0333 g/L
# Ham_medium["Ca2+[s]"] = 0.0333 * 1e3 / 110.98 # mM
Ham_medium["Ca2+[s]"] = Inf_
# Ham_medium["chloride[s]"] = 2 * 0.0333 * 1e3 / 110.98 # mM
Ham_medium["chloride[s]"] = Inf_

# Cupric Sulfate • 5H2O (CuSO4 + 5H2O)   0.0000025 g/L
# Ham_medium["Cu2+[s]"] = 0.0000025 * 1e3 / 249.685 # mM
Ham_medium["Cu2+[s]"] = Inf_
# Ham_medium["sulfate[s]"] = 0.0000025 * 1e3 / 249.685 # mM
Ham_medium["sulfate[s]"] = Inf_

# Ferrous Sulfate • 7H2O (FeSO4.7H2O)   0.000834 g/L
# Ham_medium["Fe2+[s]"] = 0.000834 * 1e3 / 151.91 # mM
Ham_medium["Fe2+[s]"] = Inf_
# Ham_medium["sulfate[s]"] += 0.000834 * 1e3 / 151.91 # mM

# Magnesium Chloride (MgCl2)   0.0576 g/L
# Ham_medium["Mg2+[s]"] = 0.0576 * 1e3 / 95.211 # mM
Ham_medium["Mg2+[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.0576 * 1e3 / 95.211 # mM

# Potassium Chloride (KCl)  0.224 g/L
# Ham_medium["K+[s]"] = 0.224 * 1e3 / 74.5513 # mM
Ham_medium["K+[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.224 * 1e3 / 74.5513 # mM

# Sodium Bicarbonate (NaHCO3)  1.176 g/L
# Ham_medium["Na+[s]"] = 1.176 * 1e3 / 84.007 # mM
Ham_medium["Na+[s]"] = Inf_
# Ham_medium["HCO3-[s]"] = 1.176 * 1e3 / 84.007 # mM
Ham_medium["HCO3-[s]"] = Inf_

# Sodium Chloride (NaCl)  7.599 g/L
# Ham_medium["Na+[s]"] += 7.599 * 1e3 / 58.44 # mM
# Ham_medium["chloride[s]"] += 7.599 * 1e3 / 58.44 # mM

# Sodium Phosphate Dibasic (anhydrous) (Na2HPO4)  0.14204 g/L
# Ham_medium["Na+[s]"] += 2 * 0.14204 * 1e3 / 141.96 # mM
# Ham_medium["Pi[s]"] = 0.14204 * 1e3 / 141.96 # mM
Ham_medium["Pi[s]"] = Inf_

# Zinc Sulfate • 7H2O (ZnSO4. 7H2O)  0.000863 g/L
# Ham_medium["zinc[s]"] = 0.000863 * 1e3 / 287.6 # mM
Ham_medium["zinc[s]"] = Inf_
# Ham_medium["sulfate[s]"] += 0.000863 * 1e3 / 287.6 # mM

# Sodium Bicarbonate  1.176 g/L
# Ham_medium["Na+[s]"] += 1.176 * 1e3 / 58.44 # mM
# Ham_medium["HCO3-[s]"] += 1.176 * 1e3 / 84.007 # mM

## ---------------------------------------------------------------
# Vitamins

# D-Biotin   0.0000073 g/L
# Ham_medium["biotin[s]"] = 0.0000073 * 1e3 / 244.31 # mM
Ham_medium["biotin[s]"] = Inf_ 

# Choline Chloride (C5H14NO.Cl) 0.01396 g/L
# Ham_medium["choline[s]"] = 0.01396 * 1e3 / 139.62 # mM
Ham_medium["choline[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.01396 * 1e3 / 139.62 # mM

# Folic Acid   0.00132 g/L
# Ham_medium["folate[s]"] = 0.00132 * 1e3 / 441.4 # mM
Ham_medium["folate[s]"] = Inf_

# myo-Inositol   0.018 g/L
# Ham_medium["inositol[s]"] = 0.018 * 1e3 / 180.16 # mM
Ham_medium["inositol[s]"] = Inf_

# Niacinamide   0.000037 g/L
# Ham_medium["nicotinamide[s]"] = 0.000037 * 1e3 / 122.12 # mM
Ham_medium["nicotinamide[s]"] = Inf_

# D-Pantothenic Acid (hemicalcium)   0.00048 g/L
# Ham_medium["pantothenate[s]"] = 0.00048 * 1e3 / 219.23 # mM
Ham_medium["pantothenate[s]"] = Inf_

# Pyridoxine • HCl   0.000062 g/L
# Ham_medium["pyridoxine[s]"] = 0.000062 * 1e3 / 169.18 # mM
Ham_medium["pyridoxine[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.000062 * 1e3 / 169.18 # mM

# Riboflavin   0.000038 g/L
# Ham_medium["riboflavin[s]"] = 0.000038 * 1e3 / 376.36 # mM
Ham_medium["riboflavin[s]"] = Inf_

# Thiamine • HCl   0.00034 g/L
# Ham_medium["thiamin[s]"] = 0.00034 * 1e3 / 265.355 # mM
Ham_medium["thiamin[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.00034 * 1e3 / 265.355 # mM

# Vitamin B12   0.00136 g/L
# Ham_medium["cobamide-coenzyme[s]"] = 0.00136 * 1e3 / 1355.365 # mM
Ham_medium["cobamide-coenzyme[s]"] = Inf_

## ---------------------------------------------------------------
# Other
# D-Glucose   1.802  g/L
Ham_medium["glucose[s]"] = 1.802 * 1e3 / 180.156 # mM

# Hypoxanthine   0.00408 g/L
# Ham_medium["hypoxanthine[s]"] = 0.00408 * 1e3 / 136.1115 # mM
Ham_medium["hypoxanthine[s]"] = Inf_

# Linoleic Acid   0.000084 g/L
# Ham_medium["linolenate[s]"] = 0.000084 * 1e3 / 280.4472 # mM
Ham_medium["linolenate[s]"] = Inf_

# Phenol Red • Na   0.0013 g/L # PH Indicator

# Putrescine • HCl   0.000161 g/L
# Ham_medium["putrescine[s]"] = 0.000161 * 1e3 / 88.15 # mM
Ham_medium["putrescine[s]"] = Inf_
# Ham_medium["chloride[s]"] += 0.000161 * 1e3 / 88.15 # mM

# Pyruvic Acid • Na   0.11 g/L # Carbon source, not included

# Thioctic Acid   0.00021 g/L
# Ham_medium["lipoic acid[s]"] = 0.00021 * 1e3 / 206.33 # mM
Ham_medium["lipoic acid[s]"] = Inf_

# Thymidine   0.00073   g/L
# Ham_medium["thymidine[s]"] = 0.00073 * 1e3 / 242.2286 # mM
Ham_medium["thymidine[s]"] = Inf_

## ---------------------------------------------------------------
# Add
# L-Glutamine   0.146 g/L
Ham_medium["glutamine[s]"] = 0.146 * 1e3 / 146.14 # mM


## ---------------------------------------------------------------
# Amino Acids

# L-Alanine   0.009 g/L
Ham_medium["alanine[s]"] = 0.009 * 1e3 / 89.09 # mM

# L-Arginine • HCl   0.211 g/L
Ham_medium["arginine[s]"] = 0.211 * 1e3 / 174.2 # mM
# Ham_medium["chloride[s]"] += 0.211 * 1e3 / 174.2 # mM

# L-Asparagine • H2O   0.01501 g/L
Ham_medium["asparagine[s]"] = 0.01501 * 1e3 / 132.12 # mM

# L-Aspartic Acid   0.0133 g/L
Ham_medium["aspartate[s]"] = 0.0133 * 1e3 / 133.11 # mM

# L-Cysteine • HCl • H2O   0.035 g/L
Ham_medium["cysteine[s]"] = 0.035 * 1e3 / 121.16 # mM
# Ham_medium["chloride[s]"] += 0.035 * 1e3 / 121.16  # mM

# L-Glutamic Acid   0.0147 g/L
Ham_medium["glutamate[s]"] = 0.0147 * 1e3 / 147.13 # mM

# L-Glutamine   0.146 g/L
Ham_medium["glutamine[s]"] += 0.146 * 1e3 / 146.14 # mM

# Glycine   0.00751 g/L
Ham_medium["glycine[s]"] = 0.00751 * 1e3 / 75.07 # mM

# L-Histidine • 3HCl • H2O   0.02096 g/L
Ham_medium["histidine[s]"] = 0.02096 * 1e3 / 155.1546 # mM
# Ham_medium["chloride[s]"] += 3 * 0.02096 * 1e3 / 155.1546 # mM

# L-Isoleucine   0.00394 g/L
Ham_medium["isoleucine[s]"] = 0.00394 * 1e3 / 131.17 # mM

# L-Leucine   0.0131 g/L
Ham_medium["leucine[s]"] = 0.0131 * 1e3 / 131.17 # mM

# L-Lysine • HCl   0.0365 g/L
Ham_medium["lysine[s]"] = 0.0365 * 1e3 / 146.19 # mM
# Ham_medium["chloride[s]"] += 0.0365 * 1e3 / 146.19 # mM

# L-Methionine   0.00448 g/L
Ham_medium["methionine[s]"] = 0.00448 * 1e3 / 149.21 # mM

# L-Phenylalanine   0.00496 g/L
Ham_medium["phenylalanine[s]"] = 0.00496 * 1e3 / 165.19 # mM

# L-Proline   0.0345 g/L
Ham_medium["proline[s]"] = 0.0345 * 1e3 / 115.13 # mM

# L-Serine   0.0105 g/L
Ham_medium["serine[s]"] = 0.0105 * 1e3 / 105.09 # mM

# L-Threonine   0.0119 g/L
Ham_medium["threonine[s]"] = 0.0119 * 1e3 / 119.1192 # mM

# L-Tryptophan   0.00204 g/L
Ham_medium["tryptophan[s]"] = 0.00204 * 1e3 / 204.23 # mM

# L-Tyrosine • 2Na • 2H2O   0.00778 g/L
Ham_medium["tyrosine[s]"] =  0.00778 * 1e3 / 181.19 # mM
# Ham_medium["Na+[s]"] += 2 * 0.00778 * 1e3 / 181.19 # mM
  
# L-Valine   0.0117 g/L
Ham_medium["valine[s]"] = 0.0117 * 1e3 / 117.151 # mM
# -

# ## Extras

Ham_medium["O2[s]"] = Inf_
Ham_medium["H2O[s]"] = Inf_
# Not found in concentration source
# TODO search real conc for this mets
# The concentration of this components was selected to ensure
# suficient growth rate compared to the experimental values
Ham_medium["alpha-tocopherol[s]"] = Inf_
Ham_medium["gamma-tocopherol[s]"] = Inf_
Ham_medium["aquacob(III)alamin[s]"] = Inf_
Ham_medium["retinoate[s]"] = Inf_

# Testing
All_components = Any["valine[s]", "putrescine[s]", "H2O[s]", "Fe2+[s]", "arginine[s]", "tryptophan[s]", 
                    "asparagine[s]", "glutamine[s]", "serine[s]", "linolenate[s]", "cysteine[s]", 
                    "HCO3-[s]", "Ca2+[s]", "gamma-tocopherol[s]", "chloride[s]", "histidine[s]", 
                    "hypoxanthine[s]", "pyridoxine[s]", "Na+[s]", "thiamin[s]", "phenylalanine[s]", 
                    "riboflavin[s]", "Cu2+[s]", "K+[s]", "lipoic acid[s]", "alanine[s]", "Mg2+[s]", 
                    "zinc[s]", "sulfate[s]", "lysine[s]", "glycine[s]", "leucine[s]", "threonine[s]", 
                    "tyrosine[s]", "aspartate[s]", "inositol[s]", "thymidine[s]", "cobamide-coenzyme[s]", 
                    "choline[s]", "proline[s]", "retinoate[s]", "glutamate[s]", "pantothenate[s]", "folate[s]", 
                    "O2[s]", "alpha-tocopherol[s]", "Pi[s]", "isoleucine[s]", "aquacob(III)alamin[s]", 
                    "glucose[s]", "nicotinamide[s]", "methionine[s]", "biotin[s]"] |> sort
@assert all(Ham_medium |> keys |> collect |> sort .== All_components)
@assert all(c >= 0.0 for (met, c) in Ham_medium)

# Print
println("\nHam's F-12")
JSON.print(Ham_medium, 4)
println()

# Saving
file = HG.HAM_MEDIUM_FILE
tagsave(file, Dict(DATA_KEY => Ham_medium))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")

