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
@quickactivate "Chemostat_Rath2017"

import MAT
using SparseArrays
using Test
using Distributions

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Human1
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;

# +
# using IJulia # Delete
# reset_IJulia_counter(th = 10_000) = (IJulia.stdio_bytes[] > th) && (IJulia.stdio_bytes[] = 0);
# -

# ---
# ## Description
# This script take the given ec_model as a template to produce a new one from the orig_model

# ---
# ## Testing process

# This take a tINIT GEM and its respective ecModel as ec template.
# So, the resulting new ecModel must be equal to the template one
function test_process()
    
    @time @testset "Generating ecModel" begin
        
        println("Testing ecModel generation")
        input_dat = wload(RepH1.COMP_FVA_HG_INPUT_FILE);
        orig_model = Ch.Utils.uncompress_model(input_dat["orig_model"]);
        ec_model = Ch.Utils.uncompress_model(input_dat["ec_model"]);

        test_gd = Human1.EcMGDC(orig_model, ec_model);
        Human1.collect_rxns_data!(test_gd; on_skip_ec = (x...) -> false);
        Human1.collect_mets_data!(test_gd);
        Human1.collect_draw_rxns!(test_gd);
        Human1.add_prot_pool_exchange!(test_gd);
        Human1.make_all_unique!(test_gd);
        
        flush(stdout)
        println(); flush(stdout);
        
        new_ec_model = Human1.build_new_model(test_gd);
        @test all(new_ec_model.S[:] |> sort .== ec_model.S[:] |> sort)
        @test all(new_ec_model.b[:] |> sort .== ec_model.b[:] |> sort)
        @test all(new_ec_model.lb[:] |> sort .== ec_model.lb[:] |> sort)
        @test all(new_ec_model.ub[:] |> sort .== ec_model.ub[:] |> sort)
        @test all(new_ec_model.mets[:] |> sort .== ec_model.mets[:] |> sort)
        @test all(new_ec_model.rxns[:] |> sort .== ec_model.rxns[:] |> sort)
        @test all(new_ec_model.subSystems[:] .|> vec .|> first |> sort .== 
                ec_model.subSystems[:] .|> vec .|> first |> sort )
    end
    return nothing
end
test_process();

# ---
# ## Build ec template

function build_ec_template()
    
    println("\nCreating ec_template")
    
    # Pick all ec Models
    data_dir = joinpath(ecG.MODEL_RAW_DATA_DIR, "ec_GEMs/models") # TODO: package this
    ec_model_files = []
    for (root, dirs, files) in walkdir(data_dir)
        for file in files
            file = joinpath(root, file)
            if basename(file) == "ecModel_batch.mat"
                push!(ec_model_files, file)
            end
        end
    end
    
    # Collect all protless rxns in the ec reference models and
    # used as allowed protless rxns
    allowed_protless_rxns = []

    # extract data from ec models
    # This model is a superset of the others
    base_model = wload(HG.BASE_MODEL_FILE)["dat"];
    base_model = Ch.Utils.uncompress_model(base_model);
    Ch.Utils.clamp_bounds!(base_model, Human1.MAX_BOUND, Human1.ZEROTH)
    Human1.try_fba(base_model, Human1.OBJ_IDER);
    
    # TODO: search protless reactions in ec_models
    ec_models = wload(ecG.EC_BASE_REFERENCE_MODELS)[DATA_KEY]
#     ec_models = ec_models[[1]] # Test
    for (i, ec_model) in ec_models |> enumerate
        @time begin
            Ch.Utils.clamp_bounds!(ec_model, Human1.MAX_BOUND, Human1.ZEROTH)
            println("\nDoing [$i/ $(length(ec_model_files))]")
            println("Base model: ", size(base_model))
            println("ec template: ", size(ec_model))
            
            gd = Human1.EcMGDC(base_model, ec_model);
            Human1.collect_rxns_data!(gd);
            Human1.collect_mets_data!(gd);
            Human1.collect_draw_rxns!(gd);
            Human1.add_prot_pool_exchange!(gd);
            Human1.make_all_unique!(gd);
            flush(stdout)
            
            base_model = Human1.build_new_model(gd);
            println("new base model:")
            Human1.try_fba(base_model, Human1.OBJ_IDER);
            
            # Collect allowed protless rxns
            protless = Human1.collect_protless(ec_model)
            protless = ec_model.rxns[protless]
            union!(allowed_protless_rxns, protless)
        end
    end
    
    ## Adding prot to protless reactions
    # protless rxns (reactions that have not a prot_ associated and
    # wasn't found as so in the ec reference models)
    intersect!(allowed_protless_rxns, base_model.rxns)
    println("\nAdd prot_* to protless rxns")
    allowed_protless_rxns = 
        [Ch.Utils.rxnindex(base_model, rxn) for rxn in allowed_protless_rxns]
    M, N = size(base_model)
    protless_rxns = trues(N)
    protless_rxns[allowed_protless_rxns] .= false # skip justified in reference models
    protless_rxns = collect(1:N)[protless_rxns]
    protless_rxns = Human1.collect_protless(base_model, protless_rxns)
    Np = length(protless_rxns)
    println("protless_rxns: ", Np , " (", round(Np/N; digits = 3) * 100, " %)" )
    
    prot_kin_stois, prot_draw_stois = Human1.prot_stois(base_model);
    prot_kin_stoi, prot_draw_stoi = (prot_kin_stois, prot_draw_stois) .|> median
    println("prot_kin_stoi: ", prot_kin_stoi)
    println("prot_draw_stoi: ", prot_draw_stoi)
    base_model = Human1.fill_protless(base_model, protless_rxns, prot_kin_stoi, prot_draw_stoi)
    
    # test (The model must only have the allowed protless rxns)
    protless_rxns = Human1.collect_protless(base_model);
    @assert all(protless_rxns |> Set == allowed_protless_rxns |> Set)
    
    println(" "^50, "\rDone: ec_template: ", size(base_model))
    Human1.try_fba(base_model, Human1.OBJ_IDER);
    return base_model
end

# ---
# ### Process brain models

# Test
if false # isfile(ecGEMs.MODEL_EC_TEMPLATE_FILE)
    # load template
    ec_template = wload(ecGEMs.MODEL_EC_TEMPLATE_FILE)[DATA_KEY];
    ec_template = Ch.Utils.uncompress_model(ec_template);
    println("\nEc template loaded: ", size(ec_template))
else
    ec_template = build_ec_template();
    to_save = Ch.Utils.compress_model(ec_template);
    file = ecGEMs.MODEL_EC_TEMPLATE_FILE
    tagsave(file, Dict(DATA_KEY => to_save));
    println(relpath(file), " created!!!, size: ", filesize(file), " bytes")
end

# load orig input models
tINIT_brain_models = wload(tIG.tINIT_RAW_BRAIN_MODELS_FILE)[DATA_KEY];

ec_brain_models = Dict()
tINIT_brain_models = tINIT_brain_models["GTEx-brain"] # Test
for (model_id, dat) in tINIT_brain_models
    
    println("\n\n ------------- Processing $model_id ------------- \n\n")

    model = dat["metnet"]
    println("Orig model size: ", size(model))
    
    gd = Human1.EcMGDC(model, ec_template);
    Human1.collect_rxns_data!(gd);
    Human1.collect_mets_data!(gd);
    Human1.collect_draw_rxns!(gd);
    Human1.add_prot_pool_exchange!(gd);
    Human1.make_all_unique!(gd);
    flush(stdout)
            
    ec_model = Human1.build_new_model(gd);
    ec_model = Ch.Utils.compress_model(ec_model)
    println("Ec model size: ", size(ec_model))
    ec_brain_models[model_id] = ec_model
    
    flush(stdout)
end

# +
# TODO: Ramake the ec inclusion workflow to process one reaction at a time, and possible 
# from several ec reference models, this will make possible to check fba at any time.
# -

# saving
file = ecG.EC_BRAIN_RAW_MODELS_FILE
tagsave(file, Dict(DATA_KEY => ec_brain_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


