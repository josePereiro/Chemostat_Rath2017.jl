import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## ----------------------------------------------------------------------
import UtilsJL: save_data, load_data, struct_to_dict, compressed_copy
import MAT

import Chemostat
import Chemostat.SteadyState: apply_bound!
import Chemostat.LP: fva_preprocess, fba
import Chemostat.Utils: MetNet, del_rxn!, compacted_model,
                        lb!, ub!, S!, rxn_mets, rxnindex,
                        rxn_str, bounds, summary, clampfields!

import Chemostat_Rath2017: RathData, MODEL1105100000
const Rd = RathData
const M = MODEL1105100000

## ----------------------------------------------------------------------
# This file is the primary input to the processing
if !isfile(M.MODEL_RAW_MAT_FILE)
    error("$(M.MODEL_RAW_MAT_FILE) not found, you must run 'scripts/MODEL1105100000/0_make_mat_file.py'")
end

## ----------------------------------------------------------------------
# LOAD MAT MODEL
# this will be the base base_model for all the processing
mat_model = MAT.matread(M.MODEL_RAW_MAT_FILE);
mat_model = mat_model[first(keys(mat_model))];
base_model = MetNet(mat_model; reshape = true);
println("Preparing Base Model: ", size(base_model))

## ----------------------------------------------------------------------
# EXCHANGES
# bkwd and fwd splatted reactions are troublemakers for EP, but they 
# are necessary to model enzymatic costs. 
# So, we leave as least as possible. We unified the exchanges (make 
# them a unique rxn), and let the in a semi-open state (intake bloked, outtake open)
# DELETING BKWD_EXCHS
bkwd_filter(x) = startswith.(x, "EX_") && endswith(x, "__bkwd")
bkwd_exchs(model) = findall(bkwd_filter , model.rxns)
exchs = bkwd_exchs(base_model)
for rxni in exchs
    del_rxn!(base_model, rxni)
end
base_model = compacted_model(base_model);
@assert isempty(bkwd_exchs(base_model))

## ----------------------------------------------------------------------
# RENAMING FWD REACTIONS
fwd_filter(x) = startswith.(x, "EX_") && endswith(x, "__fwd")
fwd_exchs(model) = findall(fwd_filter , model.rxns)
exchs = fwd_exchs(base_model)
for exch_i in exchs
    # Renaming fwd reactions (deleting the fwd part)
    base_model.rxns[exch_i] = base_model.rxns[exch_i][1:(end - length("__fwd"))]
end
@assert isempty(fwd_exchs(base_model))

## ----------------------------------------------------------------------
# Exchanges
exch_filter(x) = startswith.(x, "EX_") || startswith.(x, "DM_")
exchs = findall(exch_filter , base_model.rxns)
for exch_i in exchs
    
    # This base_model have a cost assiciated with the 
    # exchanges, we don't wants that, 
    # it is not biologically justifiable
    S!(base_model, M.COST_MET, exch_i, 0.0)
   
    # Because, this reactions are foward unbalanced (A <-> nothing)
    # positibe (+) bounds limit the outtake of the cell and
    # negative (-) bounds limit the intake.
    # Because in the Chemostat the intakes is 
    # controlled by the medium, we'll close all intakes now
    # We'll open all outtakes
    lb!(base_model, exch_i, 0.0)
    ub!(base_model, exch_i, M.ABS_MAX_BOUND)
    
end

## ----------------------------------------------------------------------
# EXCH MET MAP
# A fast way to get the exch reaction from the metabolite and viceversa

# +
exch_met_map = Dict()
for rxn in exchs
    mets = rxn_mets(base_model, rxn)
    if length(mets) == 1
        met = base_model.mets[mets[1]]
        rxn = base_model.rxns[rxn]
        
        exch_met_map[met] = rxn
        exch_met_map[rxn] = met
    end
end

# Saving
save_data(M.EXCH_MET_MAP_FILE, exch_met_map)

## ----------------------------------------------------------------------
# BASE INTAKE INFO
# The Base model will have a medium (expressed as open intake fluxes) 
# that resemble the cultivation at xi = 1 using the set up in Rath 2017 exp A. 
# So the lb of the intakes will be directly the (negative) concentration in 
# 42_MAX_UB standard medium (see Cossio's paper). Also, a few intakes, not 
# justified in the standard medium will be add based in the intakes of the 
# original model FBA analysis.
base_intake_info = Dict()
mets_map = M.load_mets_map()

# From 42_MAX_UB standard medium
for rath_met in Rd.all_mets
    # Rath ids
    conc_met_id = "c$rath_met" # Feed medium conc id
    
    # base_model id
    model_met = mets_map[rath_met]   
    exch_rxn = exch_met_map[model_met]
    
    # 42_MAX_UB standard medium
    # we take the openest version of the inteke for building the
    # base model
    conc = maximum(Rd.val(conc_met_id, Rd.ststs, 0.0)) 
    conc == 0.0 && continue
    
    lb = -M.ABS_MAX_BOUND # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end

## ----------------------------------------------------------------------
# Required open intakes from FBA analysis (a.k.a the base_model die if not open)
for rxn in ["EX_h2o_LPAREN_e_RPAREN_", 
            "EX_o2_LPAREN_e_RPAREN"]
    base_intake_info[rxn] = Dict("c" => M.ABS_MAX_CONC, "lb" => -M.ABS_MAX_BOUND) 
end

# This met is a carbon source but is required, so, 
# I restricted till garanties a credible growth rate 
# interval and variability with xi
base_intake_info["EX_adprib_LPAREN_e_RPAREN_"] = 
    Dict("c" => 1, "lb" => -M.ABS_MAX_BOUND) 

# Saving
save_data(M.BASE_INTAKE_INFO_FILE, base_intake_info)    

## ----------------------------------------------------------------------
# APPLY MEDIUM
# The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c
ξ = maximum(Rd.val(:ξ, Rd.ststs))
apply_bound!(base_model, ξ, base_intake_info);

## ----------------------------------------------------------------------
# Niklas Biomasss

# I will modified the biomass equation of MODEL1105100000 model with data
# derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
# I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
println("Applying Niklas Biomass")
niklas_biomass = M.load_niklas_biomass()
biomass_idx = rxnindex(base_model, M.BIOMASS_IDER)
base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
for (met, y) in niklas_biomass
    S!(base_model, met, biomass_idx, y)
end

## ----------------------------------------------------------------------
# ATPM demand
# Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
# pgDW/μm³ * 1e9 = pgDW/μL
# pgDW/μL * 1e6 = pgDW/L
# pgDW/L * 1e-12 = gDW/L
# atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
# # Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

# Again, we delete the cost associate with this
# reaction. It do not represent a single enzimatic reaction
S!(base_model, M.COST_MET, M.ATPM_IDER, 0.0)

# This reaction is foward defined with respect to atp
# atp + more_reacts <-> adp + more_products
# so we limit the lower bounds as the minimum atp demand 
# that the cell metabolism must fullfill
lb!(base_model, M.ATPM_IDER, atpm_flux)
println(rxn_str(base_model, M.ATPM_IDER), " ", bounds(base_model, M.ATPM_IDER))

## ----------------------------------------------------------------------
# Saving base_model
dict_ = base_model |> struct_to_dict |> compressed_copy
save_data(M.BASE_MODEL_FILE, MetNet(dict_))

## ----------------------------------------------------------------------
# FVA Preprocess
# We will reduce the bounds interval of all the reactions using the results of FVA.
# If FVA for a flux returns fva_lb == fva_lb, then the flux is blocked to lb = fva_lb, ub = fva_ub
# The method allows you to set a block eps (lb = fva_lb - eps, ub = fva_ub + eps).
# We fully blocked eps = 0, for save computations in EP.

# +
println("Base base_model summary")
summary(base_model)
println("FVA Preprocessing")
# This can take a while
# TODO use COBRA fva for speed up
fva_preprocessed_model = fva_preprocess(base_model, eps = 0, verbose = true)
##
clampfields!(fva_preprocessed_model, [:lb, :ub]; 
    zeroth = M.ZEROTH,  abs_max = M.ABS_MAX_BOUND)
println("Fva preprocessed base_model summary")

fbaout = fba(fva_preprocessed_model, M.BIOMASS_IDER)
summary(fva_preprocessed_model, fbaout)

dict_ = fva_preprocessed_model |> struct_to_dict |> compressed_copy
save_data(M.FVA_PP_BASE_MODEL_FILE, MetNet(dict_))