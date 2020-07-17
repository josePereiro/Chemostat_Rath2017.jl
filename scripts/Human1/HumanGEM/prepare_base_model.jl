# -*- coding: utf-8 -*-
# +
import DataFrames: DataFrame
import MAT
import CSV
import JSON
import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
Rd = Chemostat_Rath2017.RathData
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
# -
# This file is the primary input to the processing
if !isfile(HG.MODEL_RAW_MAT_FILE)
    error("$(HG.MODEL_RAW_MAT_FILE) not found, you must run 'make all' fisrt (see README)!!!")
end
HG.load_all_data();

# ---
# ## Description
# ---

# This script prepare the original HumanGEM model for modeling. It include some modifications to it (see comment in code) extracted from several data sources. 

# ---
# ## Prepare base base_model
# ---
# this will be the base model for all the processing

# Mat Model
mat_model = MAT.matread(HG.MODEL_RAW_MAT_FILE)["ihuman"];
Ch.Utils.reshape_mat_dict!(mat_model)
base_model = Ch.Utils.read_mat(HG.MODEL_RAW_MAT_FILE);
M, N = size(base_model)
println("\nLoaded Mat model: ", relpath(HG.MODEL_RAW_MAT_FILE))
println("Base Model: ", (M, N))

#
# ### Important Vals
# Some mets and rxns identifiers and default values that are relevant for the processing

obj_ider = "biomass_human";
bound_max_dflt = 1000
c_max_dflt = 99999 # Inf conc means that the metabolite will be never limiting the growth
atpm_ider = "HMR_3964";

# ### mets_readable_ids

# +
# we build a more human readable met ids
met_readable_ids = Dict()
for i in 1:M
    readable_id = mat_model["metNames"][i] * "[" * mat_model["comps"][mat_model["metComps"][i]] * "]";
    met_readable_ids[base_model.mets[i]] = readable_id
    met_readable_ids[readable_id] = base_model.mets[i]
end

# Saving
df = DataFrame(collect.([keys(met_readable_ids), values(met_readable_ids)]));
CSV.write(HG.BASE_READABLE_MET_IDS_FILE, df)
println("\nMet readable ids")
println("created $(relpath(HG.BASE_READABLE_MET_IDS_FILE))")
# -

# ### Exchanges
# bkwd and fwd splitted reactions are troublemakers for EP, but they are necessary to model enzymatic costs. So, we leave as least as possible. We unified the exchanges (make them a unique rxn), and let the in a semi-open state (intake bloked, outtake open)

# redefining the Chemostat exchange criteria
exch_subsys = "Any[\"Exchange/demand reactions\"]"
Ch.Utils.is_exchange(model::Ch.Utils.MetNet, ider::Ch.Utils.IDER_TYPE) = 
    exch_subsys == model.subSystems[Ch.Utils.rxnindex(model, ider)]

# I will delete all the Boundary (x) (see comps in the matmodel) metabilites, 
# leaving only the Extracellular (s) metabolites in the exchange reactions. 
# Why? they are not required
println("\nDeleting Boundary (x) metabolites ")
println("Before: ", size(base_model))
to_del = [met for met in base_model.mets if endswith(met, "x")];
base_model = Ch.Utils.del_met(base_model, to_del);
println("After: ", size(base_model))
to_del = [met for met in base_model.mets if endswith(met, "x")];
@assert isempty(to_del)

# +
# Exchanges
exchs = []
for exch_i in Ch.Utils.exchanges(base_model)
    
    # First close it, later if it is what I want, open the outtake
    Ch.Utils.lb!(base_model, exch_i, 0.0)
    Ch.Utils.ub!(base_model, exch_i, 0.0)
    
    mets = Ch.Utils.rxn_mets(base_model, exch_i)
    reacts = Ch.Utils.rxn_reacts(base_model, exch_i)
    
    length(reacts) != length(mets) != 1 && continue # I want only the forward monomoleculars
    !endswith(base_model.mets[first(reacts)], "s") && continue # I what only the exchanges 's'
    
    # Because, this reactions are forward unbalanced (A <-> nothing)
    # positibe (+) bounds limit the outtake of the cell and
    # negative (-) bounds limit the intake.
    # Because in the Chemostat the intakes is 
    # controlled by the medium, we'll close all intakes to handle with them later
    # We'll open all outtakes
    Ch.Utils.lb!(base_model, exch_i, 0.0)
    Ch.Utils.ub!(base_model, exch_i, bound_max_dflt)
    
    push!(exchs, exch_i)
    
end
println("\nExchanges: ", exchs |> length)
# -

# ### Exch Met map
# A fast way to get the exch reaction from the metabolite and viceversa

# +
exch_met_map = Dict()
for rxn in exchs
    mets = Ch.Utils.rxn_mets(base_model, rxn)
    if length(mets) == 1 # only the monomoleculars
        met = base_model.mets[mets[1]]
        rxn = base_model.rxns[rxn]
        
        exch_met_map[met] = rxn
        exch_met_map[rxn] = met
    end
end
println("\nExchanges met map: ", exch_met_map |> length)

# Saving
df = DataFrame([collect(keys(exch_met_map)), collect(values(exch_met_map))])
CSV.write(HG.EXCH_MET_MAP_FILE, df)
println("created $(relpath(HG.EXCH_MET_MAP_FILE))")
# -

# ### Base intake info
# The Base model will have a medium (expressed as open intake fluxes) that resamble the cultivation at xi = 1 using the set up in Rath 2017 exp A. So the lb of the intakes will be directly the (negative) concentration in 42_MAX_UB standard medium (see Cossio's paper). Also, a few intakes, not justified in the standard medium will be add based in the intakes of the original model FBA analysis.

# +
base_intake_info = Dict()

# From 42_MAX_UB standard medium
for rath_met in Rd.all_mets
    # Rath ids
    conc_met_id = "c$rath_met" # Feed medium conc id
    
    # base_model id
    model_met = HG.mets_map[rath_met]   
    exch_rxn = exch_met_map[model_met]
    
    # 42_MAX_UB standard medium
    # we take the openest version of the intakes for building the
    # base model
    conc = maximum(Rd.val(conc_met_id, Rd.ststs, 0.0)) 
    lb = -bound_max_dflt # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end

println("\nBase intake info: ")
JSON.print(base_intake_info, 4)
println()

# From ham's medium
for (Ham_id, conc) in HG.ham_medium
    model_met = met_readable_ids[Ham_id]
    exch_rxn = exch_met_map[model_met]
    haskey(base_intake_info, exch_rxn) && continue # Not overwrite 42_MAX_UB standard medium
    
    lb = -bound_max_dflt # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end

println("\nHams medium: ")
JSON.print(base_intake_info, 4)
println()

# # Saving
println("\nSaving")
ids = base_intake_info |> keys |> collect |> sort;
cs = [base_intake_info[id]["c"] for id in ids];
lbs = [base_intake_info[id]["lb"] for id in ids];
df = DataFrame(id = ids, c = cs, lb = lbs)
CSV.write(HG.BASE_INTAKE_INFO_FILE, df)
println("created $(relpath(HG.BASE_INTAKE_INFO_FILE))")
# -

# ### Apply medium
# The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c
ξ = 1
println("\nApplaying chemostat bound, xi: ", ξ)
Ch.SteadyState.apply_bound!(base_model, ξ, base_intake_info);


# ### Niklas Biomasss

# I will modified the biomass equation of MODEL1105100000 model with data
# derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
# I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
println("\nApplying Niklas Biomass")
JSON.print(HG.niklas_biomass, 4)
println()

biomass_idx = Ch.Utils.rxnindex(base_model, obj_ider)
base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
for (met, y) in HG.niklas_biomass
    Ch.Utils.S!(base_model, met, biomass_idx, y)
end

# ### ATPM demand

# +
# Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
# pgDW/μm³ * 1e9 = pgDW/μL
# pgDW/μL * 1e6 = pgDW/L
# pgDW/L * 1e-12 = gDW/L
# atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
# Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

# This reaction is foward defined with respect to atp
# atp + more_reacts <-> adp + more_products
# so we limit the lower bounds as the minimum atp demand 
# that the cell metabolism must fullfill
Ch.Utils.lb!(base_model, atpm_ider, atpm_flux)

# println(atpm_ider)
println("\nATP demand")
println(Ch.Utils.rxn_str(base_model, atpm_ider), " ", Ch.Utils.bounds(base_model, atpm_ider))
# -

println("\nBase Model summary")
Ch.Utils.summary(base_model)

# ---
# ## FBA Test
# ---

fbaout = Ch.LP.fba(base_model, obj_ider);
println("\nFBAout summary")
Ch.Utils.summary(base_model, fbaout)

# Saving base_model
serialize(HG.BASE_MODEL_FILE, base_model)
println("created $(relpath(HG.BASE_MODEL_FILE))!!!")


