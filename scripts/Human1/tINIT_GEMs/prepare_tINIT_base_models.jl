# -*- coding: utf-8 -*-
# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import JSON
import MAT
using SparseArrays

import Chemostat
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, HumanGEM, tINIT_GEMs, RathData
const HG = HumanGEM
const tIG = tINIT_GEMs
const Rd = RathData;
# -

# ---
# ## Description
# ---

# This script will load all the raw iINIt models and prepare then for modeling. See comments in the code for details

# ---
# ## Load brain models
# ---

tINIT_brain_models = wload(tIG.tINIT_RAW_BRAIN_MODELS_FILE)[DATA_KEY];

# ### Test Base Model

dat_id = "GTEx-brain"
println("\n ----------------- Test Processing $dat_id -----------------\n")
mat_model = tINIT_brain_models[dat_id]["mat"]
base_model = tINIT_brain_models[dat_id]["metnet"]
println("model size ", size(base_model))

#
# ### Important Vals
# Some mets and rxns identifiers and default values that are relevant for the processing

obj_ider = "biomass_human";
bound_max_dflt = 1000
c_max_dflt = 99999 # Inf conc means that the metabolite will be never limiting the growth
atpm_ider = "HMR_3964"; # From HumanGem

# ### Exchanges
# bkwd and fwd splitted reactions are troublemakers for EP, but they are necessary to model enzymatic costs. So, we leave as least as possible. We unified the exchanges (make them a unique rxn), and let the in a semi-open state (intake bloked, outtake open)        

# redefining the Chemostat exchange criteria
exch_subsys_hint = "Exchange/demand"
function is_exchange(model::Ch.Utils.MetNet, ider::Ch.Utils.IDER_TYPE)
    idx = Ch.Utils.rxnindex(model, ider)
    subsys = model.subSystems[idx]
    return occursin(exch_subsys_hint, string(subsys))
end

function delete_boundary_mets()
    # I will delete all the Boundary (x) (see comps in the matmodel) metabilites, 
    # leaving only the Extracellular (s) metabolites in the exchange reactions. 
    # Why? they are not required
    println("\nDeleting Boundary (x) metabolites ")
    println("Before: ", size(base_model))
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    global base_model = Ch.Utils.del_met(base_model, to_del);
    println("After: ", size(base_model))
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    @assert isempty(to_del)
end
delete_boundary_mets()

function define_exchanges()
    # Exchanges
    global exchs = []
    for exch in filter((ider) -> is_exchange(base_model, ider), base_model.rxns)
        exch_i = Ch.Utils.rxnindex(base_model, exch)

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
end
define_exchanges()

# ### Apply medium
# The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c


function apply_chstat_bound()
    ξ = 1.0
    println("\nApplaying chemostat bound, xi: ", ξ)
    intake_info = deepcopy(HG.base_intake_info)
    empty!(base_model.intake_info)
    for (exch, dat) in intake_info
        if !(exch in base_model.rxns)
            delete!(intake_info, exch)
        end
    end
    empty!(base_model.intake_info)
    Ch.SteadyState.apply_bound!(base_model, ξ, intake_info);
end
apply_chstat_bound();


# ### Niklas Biomasss

function apply_niklas_biomass()
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
end
apply_niklas_biomass();

# ### ATPM demand

function set_atp_demand()
    # Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
    # pgDW/μm³ * 1e9 = pgDW/μL
    # pgDW/μL * 1e6 = pgDW/L
    # pgDW/L * 1e-12 = gDW/L
    # atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
    # Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
    atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

    # From HumanGEM ATPM
    # SUMMARY (color code: warning, info, error)
    #  HMR_3964 ()
    #  lb: 0.0, ub: 1000.0
    #  (-1.0) m01371c + (-1.0) m02040c ==> (1.0) m01285c + (1.0) m02039c + (1.0) m02751c
    if !(atpm_ider in base_model.rxns)
        global base_model = Ch.Utils.add_rxn(base_model, atpm_ider; 
                mets = Dict("m01371c" => -1.0, "m02040c" => -1.0, # ->
                    "m01285c" => 1.0, "m02039c" => 1.0, "m02751c" => 1.0), 
                lb = -bound_max_dflt,
                ub = bound_max_dflt);
    end
    # This reaction is foward defined with respect to atp
    # atp + more_reacts <-> adp + more_products
    # so we limit the lower bounds as the minimum atp demand 
    # that the cell metabolism must fullfill
    Ch.Utils.lb!(base_model, atpm_ider, atpm_flux)

    # println(atpm_ider)
    println("\nATP demand")
    Ch.Utils.summary(base_model, atpm_ider)
end
set_atp_demand();

function show_model_summary()
    println("\nBase Model summary")
    Ch.Utils.summary(base_model)
end
show_model_summary();

# ---
# ## FBA Test
# ---

function run_fba_test()
    fbaout = Ch.LP.fba(base_model, obj_ider);
    println("\nFBAout summary")
    Ch.Utils.summary(base_model, fbaout)

    println("\nComparing with experiments")
    model = deepcopy(base_model)
    for stst in Rd.ststs
        println("\nStst: $stst")
        ξ = Rd.val(:ξ, stst)
        println("exp xi: $ξ")
        exp_μ = Rd.val(:μ, stst)
        println("exp growth: $exp_μ")    
        Ch.SteadyState.apply_bound!(model, ξ, Dict());
        fbaout_μ = Ch.LP.fba(model, obj_ider).obj_val;
        fbaout_μ < exp_μ && @warn("fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)")
        println("fba growth: $fbaout_μ")
    end
end
run_fba_test();

# ---
# ## Preparing  every model
# ---

# +
tINIT_base_models = Dict()
for (model_id, model) in tINIT_brain_models
    println("\n ----------------- Processing $model_id -----------------\n")
    global mat_model = tINIT_brain_models[model_id]["mat"]
    global base_model = tINIT_brain_models[model_id]["metnet"]
    println("model size ", size(base_model))

    delete_boundary_mets()
    define_exchanges()
    apply_chstat_bound()
    apply_niklas_biomass()
    set_atp_demand()
    show_model_summary()
    run_fba_test()
    
    base_model = Ch.Utils.compress_model(base_model)
    tINIT_base_models[model_id] = Dict()
    tINIT_base_models[model_id]["metnet"] = base_model
    tINIT_base_models[model_id]["mat"] = mat_model

end
# -

# Saving
file = tIG.tINIT_BASE_BRAIN_MODELS_FILE
tagsave(file, Dict(DATA_KEY => tINIT_base_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
