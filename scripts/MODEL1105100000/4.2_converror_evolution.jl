import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

@time begin
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
    const ChEP = Ch.MaxEntEP

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const M = ChR.MODEL1105100000

    import UtilsJL
    const UJL = UtilsJL

    using Plots
    import GR
    GR.inline("png")

    using Base.Threads
end

## ----------------------------------------------------------------------------
function load_model(name, stst; compressed = false)
    MINDEX = UJL.load_data(M.MODEL_INDEX_FILE; verbose = false)
    mfile = MINDEX[stst][name]
    model = deserialize(mfile)
    compressed ? model : ChU.uncompressed_model(model)
end

## ----------------------------------------------------------------------------
let
    # setup
    fig_dir = joinpath(M.MODEL_FIGURES_DATA_DIR, "4.2_err_progress")
    mkpath(fig_dir)
    
    WLOCK = ReentrantLock()

    # params
    stst = "E"
    name = :scaled_fva_pp_model
    model = load_model(name, stst)
    
    damps = [0.9, 0.95, 0.995]
    save_frec = 50
    errs_pool = Dict(damp => [] for damp in damps)
    function plot_progress()
        thid = threadid()
        for (damp, errs) in errs_pool
            try
                errs_len = length(errs)
                curr_err = last(errs)
                @info("Plotting", damp, errs_len, curr_err, now(), thid); println()
    
                p = plot(log10.(errs); title = "EP error evolution",
                    label = "", xlabel = "ep iteration", ylabel = "log(err)",
                    lw = 3
                )
                UJL.mysavefig(deepcopy(p), "err_evolution", fig_dir; name, stst, errs_len, damp)
                UJL.mysavefig(deepcopy(p), "tot_err_evolution", fig_dir; name, stst, damp)
            catch err; @warn("ERROR", damp, err); println() end
        end
    end

    # maxent_ep - err tracking
    @show save_frec damps
    @threads for damp in damps   
        thid = threadid()

        function oniter(it, epmodel)
            lock(WLOCK) do
                push!(errs_pool[damp], epmodel.stat[:max_err])
                
                errs_len = length(errs_pool[damp])
                curr_err = last(errs_pool[damp])
                @info("Doing", it, errs_len, curr_err, now(), thid); println()
            
            
                try_plotting = thid == 1 && it > 1 && rem(it, save_frec) == 0.0
                try_plotting && plot_progress()

            end
            return (false, nothing)
        end

        ChEP.maxent_ep(model; oniter, damp,
            epsconv = 1e-7, maxiter = 5000,
            verbose = false
        )
    end
end;