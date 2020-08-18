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

import MAT
using SparseArrays
using Test
using Distributions

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Rep_Human1, ecGEMs, tINIT_GEMs, HumanGEM
const RepH1 = Rep_Human1;
const ecG = ecGEMs
const tIG = tINIT_GEMs;
const HG = HumanGEM;

# +
# using IJulia # Delete
# reset_IJulia_counter(th = 10_000) = (IJulia.stdio_bytes[] > th) && (IJulia.stdio_bytes[] = 0);
# -

# ---
# ## Description
# This script take the given ec_model as a template to produce a new one from the orig_model

# ---
# ## WorkData

# +
"""
    Store all the analysis relevant data
"""
struct WorkData
    orig_model::MetNet
    ec_template::MetNet
    orig_rxns::Vector{Tuple{Int64,String}}
    orig_mets::Vector{Tuple{Int64,String}}
    ec_rxns::Vector{Tuple{Int64,String}}
    ec_mets::Vector{Tuple{Int64,String}}
    orig_rxns_left::Dict{Int64, String}
    orig_mets_left::Dict{Int64, String}
    ec_rxns_left::Dict{Int64, String}
    ec_mets_left::Dict{Int64, String}
    
    function WorkData(orig_model::MetNet, ec_template::MetNet)
        # Get cost related reactions
        # Here I'll keep all the data to be included in the ec new model
        orig_rxns = Tuple{Int64,String}[]
        orig_mets = Tuple{Int64,String}[]
        ec_rxns = Tuple{Int64,String}[]
        ec_mets = Tuple{Int64,String}[]

        # To track not included
        orig_rxns_left = Dict(i => rxn for (i, rxn) in orig_model.rxns |> enumerate) 
        orig_mets_left = Dict(i => met for (i, met) in orig_model.mets |> enumerate)
        ec_rxns_left = Dict(i => rxn for (i, rxn) in ec_template.rxns |> enumerate) 
        ec_mets_left = Dict(i => met for (i, met) in ec_template.mets |> enumerate)
        
        new(orig_model, ec_template,
            orig_rxns, orig_mets, ec_rxns, ec_mets, 
            orig_rxns_left, orig_mets_left, ec_rxns_left, ec_mets_left)
    end
end

wd_str(wd) = string("in/left ", 
                    "rxns[", 
                        "ec:",      length(wd.ec_rxns), 
                        "/",  length(wd.ec_rxns_left) , 
                        " orig:" , length(wd.orig_rxns),
                        "/" ,length(wd.orig_rxns_left),
                    "] mets[",
                        "ec:", length(wd.ec_mets), 
                        "/", length(wd.ec_mets_left), 
                        " orig:" , length(wd.orig_mets),
                        "/" , length(wd.orig_mets_left),
                    "]"
)
Base.show(io::IO, wd::WorkData) = print(io, wd_str(wd));
# -

const up_frec = 100
const prot_pool_exchange = "prot_pool_exchange"
const prot_pool = "prot_pool"

# ---
# ## Data collection functions

# IMPORTANT: this is stateful
function collect_rxns_data!(wd; verbose = true)
    verbose && println("Collecting rxns data")
    
    olen = length(wd.orig_model.rxns)
    for (oi, orxn) in wd.orig_model.rxns |> enumerate
        s1 = orxn * "No"
        s2 = orxn * "_REV"
        s3 = "arm_" * orxn
        s4 = "arm_" * orxn * "_REV"


        only_in_origin = true
        for (eci, ecrxn) in wd.ec_rxns_left
            if orxn == ecrxn || 
                startswith(ecrxn, s1) ||
                startswith(ecrxn, s2) ||
                ecrxn == s3 || 
                ecrxn == s4

                only_in_origin = false

                push!(wd.ec_rxns, (eci, ecrxn))
                delete!(wd.ec_rxns_left, eci)

                # Just printing progress
                verbose && mod(oi, up_frec) == 0 && 
                    (Core.print("[", oi, " / ", olen, "] ", wd_str(wd),
                        " ex: ", orxn, " => ", ecrxn, " "^20, "\r"); flush(stdout))

            end
        end

        if only_in_origin
            push!(wd.orig_rxns, (oi, orxn))
            delete!(wd.orig_rxns_left, oi)
        end
    end

    unique!(wd.orig_rxns)
    unique!(wd.ec_rxns)

    verbose && println(" "^50, "\r\tDone: ", wd_str(wd), " "^50)
end

function collect_mets_data!(wd; verbose = true)
    verbose && println("Collecting mets data")
    
    for (rxnsdat, metsdat, leftmet, model) in [
                                    (wd.orig_rxns, wd.orig_mets, wd.orig_mets_left, wd.orig_model),
                                    (wd.ec_rxns  , wd.ec_mets  , wd.ec_mets_left  , wd.ec_template  )]

        len = length(rxnsdat)
        for (i, (rxni, rxn)) in rxnsdat |> enumerate

            metidxs = Ch.Utils.rxn_mets(model, rxni)
            metiders = model.mets[metidxs]
            foreach(zip(metidxs, metiders)) do (metidx, metider)
                push!(metsdat, (metidx, metider))
                delete!(leftmet, metidx)
            end


            if verbose && mod(i, up_frec) == 0
                unique!(metsdat)
                 # Just printing progress
                Core.print("[", i, " / ", len, "] ", wd_str(wd), "\r")
                flush(stdout)
            end
        end
        unique!(metsdat)

    end

    verbose && println(" "^50, "\r\tDone: ", wd_str(wd), " "^50)
end

function collect_draw_rxns!(wd; verbose = true)
    # get involved prots pseudo-metabolites
    verbose && println("Collecting draw rxns prot_pool -> prot_X")
    prot_mets = filter(wd.ec_mets) do (meti, met)
        startswith(met, "prot_") && met != prot_pool
    end
    
    # get draw rxns prot_pool -> prot_X
    foreach(prot_mets) do (prot_meti, prot_met)
        draw_rxn = "draw_" * prot_met
        draw_rxni = Ch.Utils.rxnindex(wd.ec_template, draw_rxn)
        delete!(wd.ec_rxns_left, draw_rxni)
        push!(wd.ec_rxns, (draw_rxni, draw_rxn))
    end
    verbose && println("\tDone: ", wd_str(wd), " "^50)
end

function add_prot_pool_exchange!(wd; verbose = true)
    # Add prot_pool_exchange
    verbose && println("Add prot_pool_exchange")
    
    rxni = Ch.Utils.rxnindex(wd.ec_template, prot_pool_exchange)
    delete!(wd.ec_rxns_left, rxni)
    push!(wd.ec_rxns, (rxni, prot_pool_exchange))

    meti = Ch.Utils.metindex(wd.ec_template, prot_pool)
    delete!(wd.ec_mets_left, meti)
    push!(wd.ec_mets, (meti, prot_pool))
    verbose && println("\tDone: ", wd_str(wd), " "^50)
end

function make_all_unique!(wd; verbose = true)
    verbose && println("Make all unique (just for be sure)")
    
    verbose && println("\tBefore: ", wd_str(wd), " "^50)
    for (ec_dat, orig_dat) in [(wd.ec_rxns, wd.orig_rxns), 
                               (wd.ec_mets, wd.orig_mets)]
        
        unique!(ec_dat)
        unique!(orig_dat)
        ec_iders = [ider for (i, ider) in ec_dat]
        orig_dat_ = copy(orig_dat)
        empty!(orig_dat)
        for (oi, oider) in orig_dat_
            oider in ec_iders && continue
            push!(orig_dat, (oi, oider))
        end
    end
    verbose && println("\tAfter: ", wd_str(wd), " "^50)
end

# ---
# ## new ecModel build function

# TODO: build a more complite model
function build_new_model(wd)
    
    M, N = (length(wd.orig_mets) + length(wd.ec_mets), length(wd.orig_rxns) + length(wd.ec_rxns));

    S = zeros(M, N);
    b = zeros(M);
    lb = zeros(N);
    ub = zeros(N);
    rxns = Vector{String}(undef, N);
    mets = Vector{String}(undef, M);
    subSystems = Vector(undef, N);

    # This map between ider and new model idx
    rxns_newi_map = Dict(rxn => newi for (newi, (modi, rxn)) in [wd.orig_rxns; wd.ec_rxns] |> enumerate)
    mets_newi_map = Dict(met => newi for (newi, (modi, met)) in [wd.orig_mets; wd.ec_mets] |> enumerate);
    
    # Adding mets
    for (metsdat, model) in [(wd.orig_mets, wd.orig_model), 
                             (wd.ec_mets, wd.ec_template)]
        for (modmi, met) in metsdat
            newmi = mets_newi_map[met]

            # include met data
            mets[newmi] = met
            b[newmi] = Ch.Utils.b(model, modmi)
        end
    end
    
    # Adding rxns
    for (rxnsdat, model) in [(wd.orig_rxns, wd.orig_model), 
                            (wd.ec_rxns, wd.ec_template)]

        for (modri, rxn) in rxnsdat
            newri = rxns_newi_map[rxn]

            # include rxn data (including the stoichiometry)
            rxns[newri] = rxn
            subSystems[newri] = model.subSystems[modri]
            lb[newri], ub[newri] = Ch.Utils.bounds(model, modri)
            metis = Ch.Utils.rxn_mets(model, modri)
            for modmi in metis
                met = model.mets[modmi]
                newmi = mets_newi_map[met]
                S[newmi, newri] = Ch.Utils.S(model, modmi, modri)
            end
        end
    end
    
    new_ec_model = Ch.Utils.MetNet(S, b, lb, ub, rxns, mets; subSystems = subSystems);
end

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

        test_wd = WorkData(orig_model, ec_model);
        collect_rxns_data!(test_wd);
        collect_mets_data!(test_wd);
        collect_draw_rxns!(test_wd);
        add_prot_pool_exchange!(test_wd);
        make_all_unique!(test_wd);
        
        flush(stdout)
        println(); flush(stdout);
        
        new_ec_model = build_new_model(test_wd);
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

function collect_protless(model, idxs = eachindex(model.rxns))
    filter(idxs) do rxni
        metis = Ch.Utils.rxn_mets(model, rxni)
        metids = model.mets[metis]
        !any(startswith(met, "prot_") for met in metids)
    end
end

function prot_stois(ref_model)
    prot_kin_stois = []
    prot_draw_stois = []

    prot_pool_idx = Ch.Utils.metindex(ref_model, prot_pool)
    foreach(eachindex(ref_model.mets)) do meti
        met = ref_model.mets[meti]
        !startswith(met, "prot_") && return
        meti == prot_pool_idx && return

        rxnis = Ch.Utils.met_rxns(ref_model, meti)
        for rxni in rxnis
            rxn = ref_model.rxns[rxni]
            if startswith(rxn, "draw_")
                s = Ch.Utils.S(ref_model, prot_pool_idx, rxni)
                push!(prot_draw_stois, s) 
            else
                s = Ch.Utils.S(ref_model, meti, rxni)
                push!(prot_kin_stois, s)
            end
        end
    end

    return (kin_stois = prot_kin_stois, draw_stois = prot_draw_stois)
end

function fill_protless(model, protless_rxns, prot_kin_stoi, prot_draw_stoi)
    @assert prot_kin_stoi < 0.0
    @assert prot_draw_stoi < 0.0
    
    # Computing final dimentions
    nM, nN = size(model)
    println("Current dimention: ", (nM, nN))
    for rxni in protless_rxns
        rev = Ch.Utils.isrev(model, rxni)
        # One new reaction ('_REV') for each rev
        rev && (nN += 1)
        # One new prot_ met for each protless reaction
        nM += 1
        # One new 'draw_' reaction for each 'preot_' met
        nN += 1

    end
    println("Final dimention: ", (nM, nN))
    
    M, N = size(model)
    S = zeros(nM, nN);
    S[1:M, 1:N] .= model.S
    b = zeros(nM);
    b[1:M] .= model.b
    lb = zeros(nN);
    lb[1:N] .= model.lb
    ub = zeros(nN);
    ub[1:N] .= model.ub
    rxns = Vector{String}(undef, nN);
    rxns[1:N] .= model.rxns
    mets = Vector{String}(undef, nM);
    mets[1:M] .= model.mets
    subSystems = Vector(undef, nN);
    subSystems[1:N] .= model.subSystems;
    
    
    nmeti = size(model, 1) + 1
    nrxni = size(model, 2) + 1
    # This need to be stimated from template
    prot_pool_idx = Ch.Utils.metindex(model, prot_pool)

    for rxni in protless_rxns
        rev = Ch.Utils.isrev(model, rxni)
        rxn = model.rxns[rxni]

        ## foward defined
        # add prot_met
        S[nmeti, rxni] = prot_kin_stoi
        mets[nmeti] = "prot_" * rxn
        lb[rxni] = rev ? 0.0 : lb[rxni]

        # Adding draw rxn
        rxns[nrxni] = "draw_" * mets[nmeti]
        lb[nrxni], ub[nrxni] = (0.0, 1000.0)
        # (-prot_draw_stoi) prot_pool ==> (1.0) prot_met
        S[nmeti, nrxni] = 1.0
        S[prot_pool_idx, nrxni] = prot_draw_stoi
        subSystems[nrxni] = [""]
        nmeti += 1  
        nrxni += 1

        ## backward defined
        # prot_met was added in the foward reaction
        if rev
            S[:, nrxni] .= -S[:, rxni]
            lb[nrxni] = 0.0
            ub[nrxni] = abs(lb[rxni])
            rxns[nrxni] = rxns[rxni] * "_REV"
            subSystems[nrxni] = subSystems[rxni]
            nrxni += 1
        end

    end
    
    Ch.Utils.MetNet(S, b, lb, ub, rxns, mets; subSystems = subSystems);
end

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
    
    # TODO: search protless reactions in ec_models
    for (i, file) in ec_model_files |> enumerate
        @time begin
            ec_model = Ch.Utils.read_mat(file)
            println("\nDoing [$i/ $(length(ec_model_files))]")
            println("Base model: ", size(base_model))
            println("ec template: ", size(ec_model))
            
            wd = WorkData(base_model, ec_model);
            collect_rxns_data!(wd);
            collect_mets_data!(wd);
            collect_draw_rxns!(wd);
            add_prot_pool_exchange!(wd);
            make_all_unique!(wd);
            flush(stdout)
            
            base_model = build_new_model(wd);
            
            # Collect allowed protless rxns
            protless = collect_protless(ec_model)
            protless = ec_model.rxns[protless]
            union!(allowed_protless_rxns, protless)
        end
    end
    
    ## Adding prot to protless reactions
    # protless rxns (reactions that have not a prot_ associated and
    # wasn't found as so in the ec reference models)
    intersect!(allowed_protless_rxns, base_model.rxns)
    println("\nAdd prot_* to protless rexns")
    allowed_protless_rxns = 
        [Ch.Utils.rxnindex(base_model, rxn) for rxn in allowed_protless_rxns]
    M, N = size(base_model)
    protless_rxns = trues(N)
    protless_rxns[allowed_protless_rxns] .= false # skip justified in reference models
    protless_rxns = collect(1:N)[protless_rxns]
    protless_rxns = collect_protless(base_model, protless_rxns)
    Np = length(protless_rxns)
    println("protless_rxns: ", Np , " (", round(Np/N; digits = 3) * 100, " %)" )
    
    prot_kin_stois, prot_draw_stois = prot_stois(base_model);
    prot_kin_stoi, prot_draw_stoi = (prot_kin_stois, prot_draw_stois) .|> mean
    base_model = fill_protless(base_model, protless_rxns, prot_kin_stoi, prot_draw_stoi)
    
    # test (The model must only have the allowed protless rxns)
    protless_rxns = collect_protless(base_model);
    @assert all(protless_rxns |> Set == allowed_protless_rxns |> Set)
    
    println(" "^50, "\rDone: ec_template: ", size(base_model))
    return base_model
end

# ---
# ### Process brain models

if isfile(ecGEMs.MODEL_EC_TEMPLATE_FILE)
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
for (model_id, dat) in tINIT_brain_models
    
    println("\n\n ------------- Processing $model_id ------------- \n\n")

    model = dat["metnet"]
    println("Orig model size: ", size(model))
    
    wd = WorkData(model, ec_template);
    collect_rxns_data!(wd);
    collect_mets_data!(wd);
    collect_draw_rxns!(wd);
    add_prot_pool_exchange!(wd);
    make_all_unique!(wd);
    
    ec_model = build_new_model(wd);
    ec_model = Ch.Utils.compress_model(ec_model)
    println("Ec model size: ", size(ec_model))
    ec_brain_models[model_id] = ec_model
    
    flush(stdout)
end

# saving
file = ecG.EC_BRAIN_RAW_MODELS_FILE
tagsave(file, Dict(DATA_KEY => ec_brain_models))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


