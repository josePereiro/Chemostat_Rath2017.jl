# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# +
using DrWatson

using SparseArrays
import JSON

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, RathData
import Chemostat_Rath2017.Human1: Rep_Human1, ecGEMs, tINIT_GEMs, HumanGEM
const RepH1 = Rep_Human1;
const ecG = ecGEMs
const tIG = tINIT_GEMs;
const HG = HumanGEM
const Rd = RathData
# -

HG._load_all_data();

# ---
# ## Load brain models
# ---



brain_models = wload(ecG.EC_BRAIN_RAW_MODELS_FILE)[DATA_KEY];

for (model_id, model) in brain_models
    println(size(model))
end

#
# ### Important Vals
# Some mets and rxns identifiers and default values that are relevant for the processing

const obj_ider = "biomass_human";
const c_max_dflt = 9999 # Inf conc means that the metabolite will be never limiting the growth
const atpm_ider = "HMR_3964"; # From HumanGem
const prot_pool_exchange = "prot_pool_exchange";
const max_bound = 1e3
const zeroth = 1e-8;

# ### Test Base Model

dat_id = "GTEx-brain"
println("\n ----------------- Test Processing $dat_id -----------------\n")
base_model = deepcopy(brain_models[dat_id];)
base_model = Ch.Utils.uncompress_model(base_model)
Ch.Utils.clamp_bounds!(base_model, max_bound, zeroth)
println("model size ", size(base_model))

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

# +
function define_exchanges()
    # Exchanges
    global exchs = []
    global bkwd_exchs = []
    global fwd_exchs = []
    for exch in filter((ider) -> is_exchange(base_model, ider), base_model.rxns)
        exch_i = Ch.Utils.rxnindex(base_model, exch)

        # First close it, later if it is what I want, open the outtake
        Ch.Utils.lb!(base_model, exch_i, 0.0) 
        Ch.Utils.ub!(base_model, exch_i, 0.0)

        mets = Ch.Utils.rxn_mets(base_model, exch_i)

        length(mets) != 1 && continue # I want only monomoleculars
        !endswith(base_model.mets[first(mets)], "s") && continue # I what only the exchanges 's'

        # Because, this reactions are forward unbalanced (A <-> nothing)
        # positibe (+) bounds limit the outtake of the cell and
        # negative (-) bounds limit the intake.
        # Because in the Chemostat the intakes is 
        # controlled by the medium, we'll close all intakes to handle with them later
        # We'll open all outtakes
        react = Ch.Utils.rxn_reacts(base_model, exch_i)
        if isempty(react) 
            # backward defined (nothing <- A) 
            # A positive flux means intake (close)
            push!(bkwd_exchs, exch_i)
            Ch.Utils.ub!(base_model, exch_i, 0.0)
        else # forward defined (A -> nothing)
            # A positive flux means production (open)
            push!(fwd_exchs, exch_i)
            Ch.Utils.ub!(base_model, exch_i, bound_max_dflt)
        end

        push!(exchs, exch_i)
    end
    println("\nExchanges: ", exchs |> length)
    println("\tfwd_exchs (outputs): ", fwd_exchs |> length)
    println("\tbkwd_exchs (inputs): ", bkwd_exchs |> length)
    
end
define_exchanges()
# -

### Deleting bkwd exchanges
function del_bkwd_exchs()
    println("\nDeleting bkwd_exchs")
    println("Before: ", size(base_model))
    global base_model = Ch.Utils.del_rxn(base_model, bkwd_exchs);
    println("After: ", size(base_model))
    define_exchanges()
end
del_bkwd_exchs()

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
    println("Medium exchanges: (-)intake/(+)output")
    for (rxn, lb) in intake_info
        eq = Ch.Utils.rxn_str(base_model, rxn);
        bs = Ch.Utils.bounds(base_model, rxn);
        println("\t", rxn, ": ", eq, " ", bs |> collect)
    end
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

# ### Set tot prot

function set_tot_prot(model = base_model)
    Ch.Utils.bounds!(model, prot_pool_exchange, 0.0, 1000.0)
    Ch.Utils.summary(model, prot_pool_exchange)
end
set_tot_prot()

# ---
# ## FBA Test
# ---

# +
# function show_model_summary()
#     println("\nBase Model summary")
#     Ch.Utils.summary(base_model)
# end
# show_model_summary();

# +
function run_fba_test(model = base_model)
    println("\nFBAout summary")
    fbaout = Ch.LP.fba(model, obj_ider);
    fbaout_μ = fbaout.obj_val
    println("fba growth: $fbaout_μ")
#     Ch.Utils.summary(model, fbaout)

#     println("\nComparing with experiments")
#     model_ = deepcopy(model)
#     for stst in Rd.ststs
#         println("\nStst: $stst")
#         ξ = Rd.val(:ξ, stst)
#         println("exp xi: $ξ")
#         exp_μ = Rd.val(:μ, stst)
#         println("exp growth: $exp_μ")    
#         Ch.SteadyState.apply_bound!(model_, ξ, Dict());
#         fbaout_μ = Ch.LP.fba(model_, obj_ider).obj_val;
#         println("fba growth: $fbaout_μ")
#         fbaout_μ < exp_μ && @warn("fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)")
#         flush.([stdout, stderr])
#     end
end
run_fba_test();
# -

model = deepcopy(base_model);

exchs = filter((ider) -> is_exchange(model, ider), model.rxns);

foreach(model.rxns) do rxn 
    Ch.Utils.bounds!(model, rxn, -1e6, 1e6)
#     Ch.Utils.summary(model, exch)
end

set_tot_prot(model)

foreach(exchs) do rxn 
    Ch.Utils.bounds!(model, rxn, -1e6, 1e6)
#     Ch.Utils.summary(model, exch)
end

run_fba_test(model)

Ch.Utils.summary(model, "prot_pool")

function find_limiting(model::MetNet, fbaout; tol = 1e-3)
    
    limiting = []
    for (rxni, (flx, lb, ub)) in 
                    zip(fbaout.v, model.lb, model.ub) |> enumerate
        
        lb_dist = abs(lb - flx)
        ub_dist = abs(ub - flx)
        if lb_dist < tol
            push!(limiting, (rxni, :lb))
        elseif ub_dist < tol
            push!(limiting, (rxni, :ub))
        end
    end
    return limiting
end

findmax(fbaout.v)

Ch.Utils.summary(model, "m01231s")

fbaout = Ch.LP.fba(model, obj_ider);
println(fbaout.obj_val)
for (rxni, b) in find_limiting(model, fbaout)
    flx = fbaout.v[rxni]
    flx < zeroth && continue
#     if Ch.Utils.is_exchange(model, rxni) && flx < -0.0001
        rxn = model.rxns[rxni]
#         met = HG.readable_met_ids_map[HG.exch_met_map[rxn]]
        bound = (model.lb[rxni], model.ub[rxni])
        println("ex: ", rxn, 
#         " met: ", met,  
        " f: ", flx, " b: ", bound)
#     end
end



HG.readable_met_ids_map[HG.exch_met_map["HMR_9043"]]

HG.readable_met_ids_map[HG.exch_met_map["HMR_9070"]]


