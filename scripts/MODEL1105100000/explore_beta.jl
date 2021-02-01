@time begin
    import DrWatson: quickactivate
    quickactivate(@__DIR__, "Chemostat_Rath2017")

    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChEP = Ch.MaxEntEP
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
function prepare_model(ξ, stst)

    # model = UJL.load_cache(models_cache_id; verbose = false)
    model = UJL.load_data(M.FVA_PP_BASE_MODEL_FILE; verbose = false)
    model = ChU.uncompressed_model(model)
    
    isnothing(model) && error("Unable to load model!!")

    # intake info
    intake_info = M.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossio's paper, (see README)
    ChSS.apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)
    return model
end

## ----------------------------------------------------------------------------
let
    # file = "data/processed/MODEL1105100000/cache/temp_cache___4672596391543828091.jld"
    # file = "data/processed/MODEL1105100000/cache/temp_cache___11003619098011592984.jld"
    # file = "data/processed/MODEL1105100000/cache/temp_cache___14645529629190492806.jld"
    
    # best results
    # file = "data/processed/MODEL1105100000/cache/temp_cache___1659488105687387915.jld"
    # seed = deserialize(file)[:dat]
    
    stst = "E"
    ξ = Rd.val(:ξ, stst)
    model = prepare_model(ξ, stst)
    @show size(model)
    
    epout = ChEP.maxent_ep(model; 
        damp = 0.95, 
        epsconv = 1e-4, 
        alpha = Inf,
        maxvar = 1e50, minvar = 1e-50
    )
    serialize("epout_b0_$stst.jld", epout)
end
## ----------------------------------------------------------------------------
