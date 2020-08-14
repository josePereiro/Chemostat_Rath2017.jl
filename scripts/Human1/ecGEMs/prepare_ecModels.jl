# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,jl
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
using StatsBase
# import JSON

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, HumanGEM, RathData, Rep_Human1
const HG = HumanGEM
const Rd = RathData
const RepH1 = Rep_Human1;
# -

function finddups(C)
    seen = Set{eltype(C)}()
    dups = []
    for x in C
        if in(x, seen)
            push!(dups, x)
        else
            push!(seen, x)
        end
    end
    dups
end

# ---
# ## Description
# This script take the given ec_model as a template to produce a new one from the orig_model

# ---
# ## Load model

input_dat = wload(RepH1.COMP_FVA_HG_INPUT_FILE);
orig_model = Ch.Utils.uncompress_model(input_dat["orig_model"]);
ec_model = Ch.Utils.uncompress_model(input_dat["ec_model"]);

# ---
# ## Collecting ec data

up_frec = 50;

# +
# Get cost related reactions
# Here I'll keep all the data to be included in the ec new model
orig_rxns = Tuple{Int64,String}[]
orig_mets = Tuple{Int64,String}[]
ec_rxns = Tuple{Int64,String}[]
ec_mets = Tuple{Int64,String}[]

# To track not included
orig_rxns_left = Dict(i => rxn for (i, rxn) in orig_model.rxns |> enumerate) 
orig_mets_left = Dict(i => met for (i, met) in orig_model.mets |> enumerate)
ec_rxns_left = Dict(i => rxn for (i, rxn) in ec_model.rxns |> enumerate) 
ec_mets_left = Dict(i => met for (i, met) in ec_model.mets |> enumerate)

progs_str() = string("in/left ", 
                    "rxns[", 
                        "ec:",      length(ec_rxns), 
                        "/",  length(ec_rxns_left) , 
                        " orig:" , length(orig_rxns),
                        "/" ,length(orig_rxns_left),
                    "] mets[",
                        "ec:", length(ec_mets), 
                        "/", length(ec_mets_left), 
                        " orig:" , length(orig_mets),
                        "/" , length(orig_mets_left),
                    "]"
)


println("\nCollecting rxns data")
olen = length(orig_model.rxns)
for (oi, orxn) in orig_model.rxns |> sort |> enumerate
    s1 = orxn * "No"
    s2 = orxn * "_REV"
    s3 = "arm_" * orxn
    s4 = "arm_" * orxn * "_REV"
    

    only_in_origin = true
    for (eci, ecrxn) in ec_rxns_left
        if orxn == ecrxn || 
            startswith(ecrxn, s1) ||
            startswith(ecrxn, s2) ||
            ecrxn == s3 || 
            ecrxn == s4
            
            only_in_origin = false
            
            push!(ec_rxns, (eci, ecrxn))
            delete!(ec_rxns_left, eci)
            
            # Just printing progress
            mod(oi, up_frec) == 0 && 
                (Core.print("[", oi, " / ", olen, "] ", progs_str(),
                    " ex: ", orxn, " => ", ecrxn, " "^20, "\r"); flush(stdout))
            
        end
    end
    
    if only_in_origin
        push!(orig_rxns, (oi, orxn))
        delete!(orig_rxns_left, eci)
        
    end
end

orig_rxns = unique(orig_rxns)
ec_rxns = unique(ec_rxns)

println("Done: ", progs_str(), " "^50)

# +
println("\nCollecting mets data")
for (rxnsdat, metsdat, leftmet, model) in [
                                (orig_rxns, orig_mets, orig_mets_left, orig_model),
                                (ec_rxns  , ec_mets  , ec_mets_left  , ec_model  )]
    
    len = length(rxnsdat)
    for (i, (rxni, rxn)) in rxnsdat |> enumerate
        
        metidxs = Ch.Utils.rxn_mets(model, rxni)
        metiders = model.mets[metidxs]
        foreach(zip(metidxs, metiders)) do (metidx, metider)
            push!(metsdat, (metidx, metider))
            delete!(leftmet, metidx)
        end
        
       
        if mod(rxni, up_frec) == 0
            unique!(metsdat)
             # Just printing progress
            Core.print("[", i, " / ", len, "] ", progs_str(), "\r")
            flush(stdout)
        end
    end
    unique!(metsdat)
    
end

println("Done: ", progs_str(), " "^50)
# -

# get involved prots pseudo-metabolites
println("\nCollecting draw rxns prot_pool -> prot_X")
prot_mets = filter(ec_mets) do (meti, met)
    startswith(met, "prot_")
end
# get draw rxns prot_pool -> prot_X
draw_rxns = map(prot_mets) do (prot_meti, prot_met)
    draw_rxn = "draw_" * prot_met
    draw_rxni = Ch.Utils.rxnindex(ec_model, draw_rxn)
    delete!(ec_rxns_left, draw_rxni)
    push!(ec_rxns, (draw_rxni, draw_rxn))
end
println("Done: ", progs_str(), " "^50)

# Add prot_pool_exchange
println("\nAdd prot_pool_exchange")
let prot_pool_exchange = "prot_pool_exchange", prot_pool = "prot_pool"
    
    rxni = Ch.Utils.rxnindex(ec_model, prot_pool_exchange)
    delete!(ec_rxns_left, rxni)
    push!(ec_rxns, (rxni, prot_pool_exchange))
    
    meti = Ch.Utils.metindex(ec_model, prot_pool)
    delete!(ec_mets_left, meti)
    push!(ec_mets, (meti, prot_pool))
end
println("Done: ", progs_str(), " "^50)

println("\nMake all unique (just for be sure)")
orig_rxns = unique(orig_rxns)
orig_mets = unique(orig_mets)
ec_rxns = unique(ec_rxns)
ec_mets = unique(ec_mets)
println("Done: ", progs_str(), " "^50)

# ---
# ## Building new ecModel

# +
println("\nBuilding new model")
M, N = (length(ec_mets) + length(orig_mets), length(orig_rxns) + length(ec_rxns));
println("size: ", (M, N))

S = zeros(M, N);
b = zeros(M);
lb = zeros(N);
ub = zeros(N);
rxns = Vector{String}(undef, N);
mets = Vector{String}(undef, M);

# This map between ider and new model idx
rxns_newi_map = Dict(rxn => newi for (newi, (modi, rxn)) in [orig_rxns; ec_rxns] |> enumerate)
mets_newi_map = Dict(met => newi for (newi, (modi, met)) in [orig_mets; ec_mets] |> enumerate);
# -

# Adding mets
for (metsdat, model) in [(orig_mets, orig_model), 
                         (ec_mets, ec_model)]
    for (modmi, met) in metsdat
        newmi = mets_newi_map[met]

        # include met data
        mets[newmi] = met
        b[newmi] = Ch.Utils.b(model, modmi)
    end
end

# Adding rxns
for (rxndat, model) in [(orig_rxns, orig_model), 
                        (ec_rxns, ec_model)]

    for (modri, rxn) in rxndat
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
size(new_ec_model)

# ## Testing

@assert all(new_ec_model.S[:] |> sort .== ec_model.S[:] |> sort)
@assert all(new_ec_model.b[:] |> sort .== ec_model.b[:] |> sort)
@assert all(new_ec_model.lb[:] |> sort .== ec_model.lb[:] |> sort)
@assert all(new_ec_model.ub[:] |> sort .== ec_model.ub[:] |> sort)
@assert all(new_ec_model.mets[:] |> sort .== ec_model.mets[:] |> sort)
@assert all(new_ec_model.rxns[:] |> sort .== ec_model.rxns[:] |> sort)
