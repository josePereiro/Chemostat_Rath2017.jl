using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin
    import DataFrames: DataFrame
    import MAT
    import CSV

    import Chemostat
    const Ch = Chemostat
    import Chemostat_Rath2017: Human1, RathData
    const Rd = RathData
    const H1 = Human1
    const HG = H1.HumanGEM
end

## ---------------------------------------------------------------------
## Description

# This script prepare the original HumanGEM model for modeling. It include some modifications to it (see comment in code) extracted from several data sources. 

## ---------------------------------------------------------------------
# Prepare base base_model
# this will be the base model for all the processing

# Mat Model
@time begin
    mat_model = HG.load_humangem_mat_model()
    base_model = HG.load_humangem_raw_model(mat_model)
    base_model = Ch.Utils.uncompressed_model(base_model)
    M, N = size(base_model)
    println("\nLoaded Mat model")
    println("Base Model: ", (M, N))
end

## ---------------------------------------------------------------------
# mets_readable_ids
# we build a more human readable met ids
println("\nMet readable ids")
met_readable_ids = Dict()
for i in 1:M
    readable_id = mat_model[:metNames][i] * "[" * mat_model[:comps][mat_model[:metComps][i]] * "]";
    met_readable_ids[base_model.mets[i]] = readable_id
    met_readable_ids[readable_id] = base_model.mets[i]
end

# Saving
sdat(HG, met_readable_ids, 
    "met_readable_ids", ".jls"; 
    verbose = true
)

## ---------------------------------------------------------------------
# Exchanges
# bkwd and fwd splitted reactions are troublemakers for EP, but they are necessary to model enzymatic costs.
# So, we leave as least as possible. We unified the exchanges (make them a unique rxn),
# and let them in a semi-open state (intake bloked, outtake open)

# I will delete all the Boundary (x) (see comps in the matmodel) metabilites, 
# leaving only the Extracellular (s) metabolites in the exchange reactions. 
# Why? they are not required
println("\nDeleting Boundary (x) metabolites ")
to_del = [met for met in base_model.mets if endswith(met, "x")];
println("todel: ", length(to_del))
Ch.Utils.del_met!(base_model, to_del);

exch_subsys_hint = "Exchange/demand"
# redefining the Chemostat exchange criteria
exch_subsys_hint = "Exchange/demand"
function is_exchange(model::Ch.Utils.MetNet, ider::Ch.Utils.IDER_TYPE)
    idx = Ch.Utils.rxnindex(model, ider)
    subsys = model.subSystems[idx]
    return occursin(exch_subsys_hint, string(subsys))
end

## ---------------------------------------------------------------------
# Exchanges
exchs = []
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
    Ch.Utils.ub!(base_model, exch_i, HG.ABS_MAX_BOUND)
    
    push!(exchs, exch_i)
    
end
println("\nExchanges: ", exchs |> length)

## ---------------------------------------------------------------------
# Exch Met map
# A fast way to get the exch reaction from the metabolite and viceversa

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
sdat(HG, exch_met_map, 
    "exch_met_map", ".jls"; 
    verbose = true
)

## ---------------------------------------------------------------------
# Base intake info
# The Base model will have a medium (expressed as open intake fluxes) that resamble the cultivation
# at xi = 1 using the set up in Rath 2017 exp A. So the lb of the intakes will be directly the (negative)
# concentration in 42_MAX_UB standard medium (see Cossio's paper). Also, a few intakes, not justified in
# the standard medium will be add based in the intakes of the original model FBA analysis.

base_intake_info = Dict()
mets_map = HG.load_mets_map()

# From 42_MAX_UB standard medium
for rath_met in Rd.ALL_METS
    
    # base_model id
    model_met = mets_map[rath_met]   
    exch_rxn = exch_met_map[model_met]
    
    # 42_MAX_UB standard medium
    # we take the openest version of the intakes for building the
    # base model
    conc = maximum(Rd.cval(rath_met, Rd.STSTS, 0.0)) 
    lb = -HG.ABS_MAX_BOUND # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end

# From ham's medium see Ham_medium.jl
for (Ham_id, conc) in HG.load_Ham_medium()
    model_met = met_readable_ids[Ham_id]
    exch_rxn = exch_met_map[model_met]
    haskey(base_intake_info, exch_rxn) && continue # Not overwrite 42_MAX_UB standard medium
    
    lb = -HG.ABS_MAX_BOUND # intake bound
    base_intake_info[exch_rxn] = Dict("c" => conc, "lb" => lb) 
end

# Saving
sdat(HG, base_intake_info, 
    "base_intake_info", ".jls"; 
    verbose = true
)

## ---------------------------------------------------------------------
# Apply medium
# The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c
ξ = 1
println("\nApplaying chemostat bound, xi: ", ξ)
Ch.SteadyState.apply_bound!(base_model, ξ, base_intake_info);


## ---------------------------------------------------------------------
# Niklas Biomasss

# I will modified the biomass equation of MODEL1105100000 model with data
# derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
# I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
println("\nApplying Niklas Biomass")

biomass_idx = Ch.Utils.rxnindex(base_model, HG.HUMAN_BIOMASS_IDER)
base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
for (met, y) in HG.load_Niklas_biomass()
    Ch.Utils.S!(base_model, met, biomass_idx, y)
end

## ---------------------------------------------------------------------
# ATPM demand

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
Ch.Utils.lb!(base_model, HG.HUMAN_ATPM_IDER, atpm_flux)

println("\nATP demand")
println(Ch.Utils.rxn_str(base_model, HG.HUMAN_ATPM_IDER), " ", Ch.Utils.bounds(base_model, HG.HUMAN_ATPM_IDER))

## ---------------------------------------------------------------------
# summary
println("\nBase Model summary")
Ch.Utils.summary(base_model)

## ---------------------------------------------------------------------
#  FBA Test

fbaout = Ch.LP.fba(base_model, HG.HUMAN_BIOMASS_IDER);
println("\nFBAout summary")
Ch.Utils.summary(base_model, fbaout)

## ---------------------------------------------------------------------
println("\nComparing with experiments")
for stst in Rd.STSTS
    model = deepcopy(base_model)
    println("\nStst: $stst")
    ξ = Rd.val(:ξ, stst)
    println("exp xi: $ξ")
    exp_μ = Rd.val(:μ, stst)
    println("exp growth: $exp_μ")    
    Ch.SteadyState.apply_bound!(model, ξ, base_intake_info);
    fbaout_μ = Ch.LP.fba(model, HG.HUMAN_BIOMASS_IDER).obj_val;
    fbaout_μ < exp_μ && @warn("fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)")
    println("fba growth: $(fbaout_μ)")
    
end

## ---------------------------------------------------------------------
# Saving base_model
# Saving
sdat(HG, Ch.Utils.compressed_model(base_model), 
    "HumanGEM_base_model", ".jls"; 
    verbose = true
)