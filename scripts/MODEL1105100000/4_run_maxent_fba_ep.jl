import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w"
        help = "number of workers to use"
        default = "1"
    "--clear-at-init"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--clear-at-finish"
        help = "clear cache at the end"   
        action = :store_true
end

if isinteractive()
    # Dev values
    wcount = 0
    init_clear_flag = false
    finish_clear_flag = false
else
    parsed_args = parse_args(set)
    wcount = parse(Int, parsed_args["w"])
    init_clear_flag = parsed_args["clear-at-init"]
    finish_clear_flag = parsed_args["clear-at-finish"]
end

## ----------------------------------------------------------------------------
using Distributed

NO_WORKERS = min(length(Sys.cpu_info()), wcount)
length(workers()) < NO_WORKERS && 
    addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

## Loading everywhere
@everywhere begin

    import DrWatson: quickactivate
    quickactivate(@__DIR__, "Chemostat_Rath2017")

    using Distributed
    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSU = Ch.SimulationUtils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    
    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const M = ChR.MODEL1105100000

    import UtilsJL
    const UJL = UtilsJL
    
    UJL.set_cache_dir(M.MODEL_CACHE_DATA_DIR)
end

## ----------------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    UJL.tagprintln_inmw("CLEARING CACHE ")
    UJL.delete_temp_caches()
    UJL.println_inmw("\n")
end

## ----------------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "MAXENT_FBA_EP_v1"

## ----------------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function load_model(name, stst; compressed = false)
    MINDEX = UJL.load_data(M.MODEL_INDEX_FILE)
    mfile = MINDEX[stst][name]
    model = deserialize(mfile)
    compressed ? model : ChU.uncompressed_model(model)
end

## ----------------------------------------------------------------------------
# RES IDS
# Collect all the computed results ids for bundling
const chnl = RemoteChannel() do
    Channel{Any}(10)
end
const res_ids = []
const collector = @async while true
    id = take!(chnl)
    push!(res_ids, id)
end

## ----------------------------------------------------------------------------
# BETA EXPLORATION SIMULATION
# Any of the loops can be parallelized by just 
# changing one of the 'map' functions

# ststs = Rd.ststs
ststs = ["E"] # Test
pmap(ststs) do stst
    
    ## ----------------------------------------------------------------------------
    ## SIMULATION PARAMS
    βs = [0.0; UJL.logspace(1, 9, 1000)] 

    # This determine how often EP results will be cached
    epochlen = 30

    epconv_kwargs = Dict()
    # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:maxiter] = Int(1e3) 
    # The error threshold of convergence (Big at first)
    epconv_kwargs[:epsconv] = 1e-4
    epconv_kwargs[:damp] = 0.98
    epconv_kwargs[:maxvar] = 1e50
    epconv_kwargs[:minvar] = 1e-50

    epmodel_kwargs = Dict()
    epmodel_kwargs[:alpha] = Inf

    # STOP CONDITION
    comp_model = load_model(:scaled_fva_pp_model, stst; compressed = true)
    exp_μ = Rd.val(:μ, stst)
    μ_convth = 1e-2
    

    ## ----------------------------------------------------------------------------
    # BETAS
    sim_id = (stst, sim_global_id)
    get_model() = load_model(:scaled_fva_pp_model, stst)

    function on_betaiter(epout)
        isnothing(epout) && return (false, epout)
        
        curr_μ = ChU.av(comp_model, epout, M.BIOMASS_IDER)
        isnan(curr_μ) && error("curr_μ = NaN")
        
        err = abs(curr_μ - exp_μ)/ exp_μ
        break_ = err < μ_convth || curr_μ > exp_μ

        UJL.tagprintln_inmw("BEFORE EPOCH", 
            "\nsim id:                      ", sim_id, 
            "\ncurr_μ:                      ", curr_μ, 
            "\nexp_μ:                       ", exp_μ, 
            "\ngrowth err:                  ", err,
            "\nbreak:                         ", break_,
            "\n"
        )
        
        return break_
    end

    sim_dat = ChSU.cached_simulation(;
        sim_id, epochlen, 
        verbose = true,
        get_model, on_betaiter,
        objider = M.BIOMASS_IDER, 
        costider = M.COST_IDER,
        beta_info = [(M.BIOMASS_IDER, βs)],
        clear_cache = false, use_seed = true,
        epmodel_kwargs, epconv_kwargs
    )

    serialize("test_sim_dat.jld", sim_dat)
    
    # ## SAVING DATA
    # res_id = (:RESULT, sim_id)
    # UJL.save_cache(res_id, (stst, βs, sim_dat); 
    #     headline = "CATCHING RESULTS\n"
    # )
    
    # ## PASSING ID TO MASTER
    # put!(chnl, res_id)
    
    # GC.gc()
    return nothing
end # map(Rd.ststs) do stst

# ## ----------------------------------------------------------------------------
# # COLLECTING RESULTS
# UJL.tagprintln_inmw("COLLECTING RESULTS ")
# sleep(1) # wait for collector to get all ids
# bundles = Dict()
# for id in res_ids

#     stst, ξ, βs, model, dat = UJL.load_cache(id; verbose = false)
    
#     # Bundle
#     bundle = get!(bundles, stst, ChU.ChstatBundle())

#     bundle[ξ, :net] = model
#     bundle[ξ, :ChLP.fba] = dat[:ChLP.fba]

#     for (βi, β) in βs |> enumerate
#         bundle[ξ, β, :ep] = dat[(:ep, βi)]
#     end

# end

# ## ----------------------------------------------------------------------------
# # SAVING
# UJL.tagprintln_inmw("SAVING RESULTS ")
# UJL.save_data(M.MAXENT_FBA_EB_BUNDLES_FILE, bundles)

# ## ----------------------------------------------------------------------------
# ## CLEAR CACHE (WARNING)
# if finish_clear_flag
#     UJL.tagprintln_inmw("CLEARING CACHE ")
#     UJL.delete_temp_caches()
# end