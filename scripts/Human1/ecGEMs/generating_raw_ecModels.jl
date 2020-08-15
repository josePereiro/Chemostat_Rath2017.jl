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

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Rep_Human1, ecGEMs, tINIT_GEMs
const RepH1 = Rep_Human1;
const ecG = ecGEMs
const tIG = tINIT_GEMs;
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
    up_frec::Int
    
    function WorkData(orig_model::MetNet, ec_template::MetNet; up_frec = 50)
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
            orig_rxns_left, orig_mets_left, ec_rxns_left, ec_mets_left, 
            up_frec)
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
                verbose && mod(oi, wd.up_frec) == 0 && 
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
    
    all_mets = [] # To enforce uniqueness
    
    for (rxnsdat, metsdat, leftmet, model) in [
                                    (wd.orig_rxns, wd.orig_mets, wd.orig_mets_left, wd.orig_model),
                                    (wd.ec_rxns  , wd.ec_mets  , wd.ec_mets_left  , wd.ec_template  )]

        len = length(rxnsdat)
        for (i, (rxni, rxn)) in rxnsdat |> enumerate

            metidxs = Ch.Utils.rxn_mets(model, rxni)
            metiders = model.mets[metidxs]
            foreach(zip(metidxs, metiders)) do (metidx, metider)
                metider in all_mets && return
                push!(metsdat, (metidx, metider))
                push!(all_mets, metider)
                delete!(leftmet, metidx)
            end


            if verbose && mod(rxni, wd.up_frec) == 0
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
        startswith(met, "prot_")
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
    let prot_pool_exchange = "prot_pool_exchange", prot_pool = "prot_pool"

        rxni = Ch.Utils.rxnindex(wd.ec_template, prot_pool_exchange)
        delete!(wd.ec_rxns_left, rxni)
        push!(wd.ec_rxns, (rxni, prot_pool_exchange))

        meti = Ch.Utils.metindex(wd.ec_template, prot_pool)
        delete!(wd.ec_mets_left, meti)
        push!(wd.ec_mets, (meti, prot_pool))
    end
    verbose && println("\tDone: ", wd_str(wd), " "^50)
end

function make_all_unique!(wd; verbose = true)
    verbose && println("Make all unique (just for be sure)")
    unique!(wd.orig_rxns)
    unique!(wd.orig_mets)
    unique!(wd.ec_rxns)
    unique!(wd.ec_mets)
    verbose && println("\tDone: ", wd_str(wd), " "^50)
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
            lb[newri], ub[newri] = Ch.Utils.bounds(model, modri)
            metis = Ch.Utils.rxn_mets(model, modri)
            for modmi in metis
                met = model.mets[modmi]
                newmi = mets_newi_map[met]
                S[newmi, newri] = Ch.Utils.S(model, modmi, modri)
            end
        end
    end
    
    new_ec_model = Ch.Utils.MetNet(S, b, lb, ub, rxns, mets);
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
        println()
        
        new_ec_model = build_new_model(test_wd);
        @test all(new_ec_model.S[:] |> sort .== ec_model.S[:] |> sort)
        @test all(new_ec_model.b[:] |> sort .== ec_model.b[:] |> sort)
        @test all(new_ec_model.lb[:] |> sort .== ec_model.lb[:] |> sort)
        @test all(new_ec_model.ub[:] |> sort .== ec_model.ub[:] |> sort)
        @test all(new_ec_model.mets[:] |> sort .== ec_model.mets[:] |> sort)
        @test all(new_ec_model.rxns[:] |> sort .== ec_model.rxns[:] |> sort)
    end
    return nothing
end
test_process()

# ---
# ## generate brain ecModels

# +
# load template
ec_template = wload(ecGEMs.MODEL_EC_TEMPLATE_FILE)[DATA_KEY];
ec_template = Ch.Utils.uncompress_model(ec_template);

# load orig input models
tINIT_brain_models = wload(tIG.tINIT_BRAIN_MODELS_FILE)[DATA_KEY];
# -

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


