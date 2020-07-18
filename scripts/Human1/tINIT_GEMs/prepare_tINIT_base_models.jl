# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
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
import JSON
import DataFrames: DataFrame
import MAT
import CSV
import Serialization: serialize, deserialize

import Chemostat
Ch = Chemostat
import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM
HG.load_all_data()

iTG = Chemostat_Rath2017.tINIT_GEMs

Rd = Chemostat_Rath2017.RathData
# Rd.load_all_data()
# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();
# -

# ---
# ## Description
# ---

# This script will load all the raw iINIt models and prepare then for modeling. See comments in the code for details

# ---
# ## Mat iINIT output files
# ---

healty_models_id = ["brain", "GBM NT"];
cancer_models_id = ["brain-cancer", "GBM TP", "GBM TR", "LGG TP", "LGG TR"];

raw_data_dir = "/Users/Pereiro/University/Research/Metabolism/MaxEntEP2020/"*
    "Working_Version/Chemostat_Rath2017/data/raw"
run_tINIT_outputs = joinpath(raw_data_dir, 
        "Human1_Publication_Data_Scripts/tINIT_GEMs/run_tINIT_outputs");

mat_files = Dict()
mat_files["GTEx"] = joinpath(run_tINIT_outputs, "GTEx/tINIT_GTEx_outputs.mat");
mat_files["Hart2015"] = joinpath(run_tINIT_outputs, "Hart2015/tINIT_Hart2015_HumanGEM_outputs.mat");
mat_files["TCGA"] = joinpath(run_tINIT_outputs, "TCGA/tINIT_TCGA_outputs.mat")
mat_files["DepMap_1_1"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_1_1.mat")
mat_files["DepMap_1_2"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_1_2.mat")
mat_files["DepMap_2_1"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_2_1.mat")
mat_files["DepMap_2_2"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_2_2.mat")
mat_files["DepMap_3_1"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_3_1.mat")
mat_files["DepMap_3_2"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_3_2.mat")
mat_files["DepMap_4_1"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_4_1.mat")
mat_files["DepMap_4_2"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_4_2.mat")
mat_files["DepMap_5_1"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_5_1.mat")
mat_files["DepMap_5_2"] = joinpath(run_tINIT_outputs, "DepMap_HumanGEM/tINIT_DepMap_HumanGEM_outputs_5_2.mat")
dat_ids = mat_files |> keys |> collect;

function get_brain_models(dat_key, dat)
    brain_models = Dict()
    mat_models = dat["model"]
    ids = get(dat, "id", get(dat, "tissues", nothing))
    isnothing(ids) && error("ids not found!!!")
    for (id, mat_model) in zip(ids, mat_models)
        if id in healty_models_id || id in cancer_models_id
            model_key = "$dat_key-$id"
            model_dict = get!(brain_models, model_key, Dict())
            mat_ = get!(model_dict, "mat", [])
            metnet_ = get!(model_dict, "metnet", [])
            push!(mat_, mat_model)
            push!(metnet_, Ch.Utils.MetNet(mat_model))
        end
    end
    return brain_models
end

# ### Test Base Model

dat_id = "GTEx"
@time mat_dat = MAT.matread(mat_files[dat_id])["INIT_output"]
brain_models = get_brain_models(dat_id, mat_dat)
first_key = brain_models |> keys |> collect |> first
println("\n ----------------- Test Processing $first_key -----------------\n")
mat_model = brain_models[first_key]["mat"][1];
base_model = brain_models[first_key]["metnet"][1];
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
    Ch.SteadyState.apply_bound!(base_model, ξ, HG.base_intake_info);
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
        Ch.SteadyState.apply_bound!(model, ξ, HG.base_intake_info);
        fbaout = Ch.LP.fba(model, obj_ider);
        println("fba growth: $(fbaout.obj_val)")
        
    end
end
run_fba_test();



# ---
# ## Iterate for every model
# ---

# TODO: parallelize this
for (dat_id, file) in mat_files
    @time mat_dat = MAT.matread(file)["INIT_output"]
    println("\n\n ----------------- $dat_id loaded -----------------\n")
    println("content: ", mat_dat |> keys |> collect, "\n")
    flush(stdout)
    
    brain_models = get_brain_models(dat_id, mat_dat)
    for (id, dat) in brain_models
        for (i,(metnet, mat)) in zip(dat["metnet"], dat["mat"]) |> enumerate
            model_id = "$id-$i"
            println("\n ----------------- Processing $model_id -----------------\n")
            global base_model = metnet
            println("model size ", size(base_model))
            
            delete_boundary_mets()
            define_exchanges()
            apply_chstat_bound()
            apply_niklas_biomass()
            set_atp_demand()
            show_model_summary()
            run_fba_test()
            
            # Saving
            file = "$model_id.jls"
            serialize(file, (model_id, base_model, mat))
            # TODO: package this
            println(relpath(file), " created")
        end
    end
    
    break; #Test
end;


