# -*- coding: utf-8 -*-
# +
import DataFrames: DataFrame
import MAT
import CSV
import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
M = Chemostat_Rath2017.MODEL1105100000
Rd = Chemostat_Rath2017.RathData

# -
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# This file is the primary input to the processing
if !isfile(M.MODEL_RAW_MAT_FILE)
    error("$(M.MODEL_RAW_MAT_FILE) not found, you must run 'make all' fisrt (see README)!!!")
end


# ---
# ## Prepare base base_model
# ---
# this will be the base base_model for all the processing

mat_model = MAT.matread(M.MODEL_RAW_MAT_FILE);
mat_model = mat_model[first(keys(mat_model))];
base_model = Ch.Utils.MetNet(mat_model);
println("Preparing Base Model: ", size(base_model))

#
# ### Important Vals
# Some mets and rxns identifiers oand default values that are relevant for the processing

obj_ider = "BIOMASS";
cost_met = "enzyme_c";
bound_max_dflt = 1000
c_max_dflt = 1e20 # Inf conc means that the metabolite will be never limiting the growht
atpm_ider = "ATPmaint_LSQBKT_c_RSQBKT_";

# ### Exchanges
# bkwd and fwd splitted reactions are troublemakers for EP, but they are necessary to base_model metabolic costs. So, we leave as least as possible. We unified the exchanges (make them a unique rxn), and let the in a semi-open state (intake bloked, outtake open)

# +
# Deleting bkwd_exchs
bkwd_filter(x) = startswith.(x, "EX_") && endswith(x, "__bkwd")
base_model = Ch.Utils.del_rxn(base_model, findall(bkwd_filter , base_model.rxns));
@assert isempty(findall(bkwd_filter , base_model.rxns))

# Renaming fwd reactions
fwd_filter(x) = startswith.(x, "EX_") && endswith(x, "__fwd")
fwd_exchs = findall(fwd_filter , base_model.rxns)
for exch_i in fwd_exchs
    # Renaming fwrd reactions (deleting the fwd part)
    base_model.rxns[exch_i] = base_model.rxns[exch_i][1:(end - length("__fwd"))]
end
@assert isempty(findall(fwd_filter , base_model.rxns))

# Exchanges
exch_filter(x) = startswith.(x, "EX_") || startswith.(x, "DM_")
exchs = findall(exch_filter , base_model.rxns)
for exch_i in exchs
    
    # This base_model have a cost assiciated with the 
    # exchanges, we don't wants that, 
    # it is not biologically justifiable
    Ch.Utils.S!(base_model, cost_met, exch_i, 0.0)
   
    # Because, this reactions are foward unbalanced (A <-> nothing)
    # positibe (+) bounds limit the outtake of the cell and
    # negative (-) bounds limit the intake.
    # Because in the Chemostat the intake is 
    # controlled by the medium, we'll close now
    # We'll open all outtakes
    Ch.Utils.lb!(base_model, exch_i, 0.0)
    Ch.Utils.ub!(base_model, exch_i, bound_max_dflt)
    
end
# -

# ### Exch Rxn map
# A fast way to get the exch reaction from the metabolite and viceversa

# +
exch_met_map = Dict()
for rxn in exchs
    mets = Ch.Utils.rxn_mets(base_model, rxn)
    if length(mets) == 1
        met = base_model.mets[mets[1]]
        rxn = base_model.rxns[rxn]
        
        exch_met_map[met] = rxn
        haskey(exch_met_map, rxn) && error("Duplicated key")
        exch_met_map[rxn] = met
    end
end

# Saving
df = DataFrame([collect(keys(exch_met_map)), collect(values(exch_met_map))])
CSV.write(M.EXCH_MET_MAP_FILE, df)
println("created $(relpath(M.EXCH_MET_MAP_FILE))")
# -

# ### Base intake info
# The Base base_model will have a medium (expressed as open intake fluxes) that resamble the cultivation at xi = 1 using the set up in Rath 2017 exp A. So the lb of the intakes will be directly the (the negative) concentration in 42_MAX_UB standard medium. Also, a few intakes, not justified in the standart medium will be add based in the intakes of the original base_model FBA analysis.

# +
base_intake_info = Dict()

# From 42_MAX_UB standard medium
for rath_met in Rd.all_mets
    # Rath ids
    conc_met_id = "c$rath_met" # Feed medium conc id
    
    # base_model id
    model_met = M.mets_map[rath_met]   
    exch_rxn = exch_met_map[model_met]
    
    # 42_MAX_UB standard medium
    # we take the openest version of the inteke for building the
    # base model
    conc = maximum(Rd.val(conc_met_id, Rd.ststs, 0.0)) 
    conc == 0.0 && continue
    
    lb = -bound_max_dflt # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end


# Required open intakes from FBA analysis (a.k.a the base_model die if not open)
for rxn in ["EX_h2o_LPAREN_e_RPAREN_", 
            "EX_o2_LPAREN_e_RPAREN"]
    base_intake_info[rxn] = Dict("c" => c_max_dflt, "lb" => -bound_max_dflt) 
end

# This met is a carbon source but is required, so, 
# I restricted till garanties a credible growth rate 
# interval and variability with xi
base_intake_info["EX_adprib_LPAREN_e_RPAREN_"] = 
    Dict("c" => 1, "lb" => -bound_max_dflt) 


# Saving
ids = base_intake_info |> keys |> collect |> sort;
cs = [base_intake_info[id]["c"] for id in ids];
lbs = [base_intake_info[id]["lb"] for id in ids];
df = DataFrame(id = ids, c = cs, lb = lbs)
CSV.write(M.BASE_INTAKE_INFO_FILE, df)
println("created $(relpath(M.BASE_INTAKE_INFO_FILE))")
# -

# ### Apply medium
# The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c

ξ = 1
Ch.SteadyState.apply_bound!(base_model, ξ, base_intake_info);


# ### Niklas Biomasss

# +
# I will modified the biomass equation of MODEL1105100000 model with data
# derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
# I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
println("Applying Niklas Biomass")
biomass_idx = Ch.Utils.rxnindex(base_model, obj_ider)
base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
for (met, y) in M.niklas_biomass
    met_idx = Ch.Utils.metindex(base_model, met)
    base_model.S[met_idx, biomass_idx] = y
end


# ### ATPM demand

# +
# Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
# pgDW/μm³ * 1e9 = pgDW/μL
# pgDW/μL * 1e6 = pgDW/L
# pgDW/L * 1e-12 = gDW/L
# atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
# # Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

# Again, we delete the cost associate with this
# reaction. It do not represent a single enzimatic reaction
Ch.Utils.S!(base_model, cost_met, atpm_ider, 0.0)

# This reaction is foward defined with respect to atp
# atp + more_reacts <-> adp + more_products
# so we limit the lower bounds as the minimum atp demend 
# that the cell metabolism must fullfill
Ch.Utils.lb!(base_model, atpm_ider, atpm_flux)

# println(atpm_ider)
# println(Ch.Utils.rxn_str(base_model, atpm_ider), " ", Ch.Utils.bounds(base_model, atpm_ider))
# -

# Saving base_model
serialize(M.BASE_MODEL_FILE, base_model)
println("created $(relpath(M.BASE_MODEL_FILE))!!!")

# ---
# ## FVA Preprocess
# ---
# We will reduce the bounds interval of all the reactions using the results of FVA.
# If FVA for a flux returns fva_lb == fva_lb, then the flux is blocked to lb = fva_lb, ub = fva_ub
# The method allows you to set a block eps (lb = fva_lb - eps, ub = fva_ub + eps).
# We fully blocked eps = 0, for save computations in EP.

# +
println("Base base_model summary")
Ch.Utils.summary(base_model)

println("FVA Preprocessing")
# This can take a while
# TODO use COBRA fva for speed up
fva_preprocessed_model = Ch.Utils.preprocess(base_model, eps = 0, verbose = true)

println("Fva preprocessed base_model summary")
Ch.Utils.summary(fva_preprocessed_model)
serialize(M.FVA_PP_BASE_MODEL_FILE, fva_preprocessed_model)
println("$(relpath(M.FVA_PP_BASE_MODEL_FILE)) created")
# -