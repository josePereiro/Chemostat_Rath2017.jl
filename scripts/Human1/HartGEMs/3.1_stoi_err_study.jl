using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using ProgressMeter
    using Base.Threads

    import Chemostat
    import Chemostat.MetNets
    import Chemostat.MetLP
    import Chemostat.MetEP
    const Ch = Chemostat

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
    const AG = H1.HartGEMs

    using Plots
    using Statistics
end

## ---------------------------------------------------------------------
function _find_epout_files(modelid, tissue)
    modelid = string(modelid)
    tissue = string(tissue)
    files = Dict()
    for file in readdir(datdir(AG); join = true)
        heads, params, ext = parse_dfname(file)
        !("maxent_ep_beta0" in heads) && continue
        !("epout" in heads) && continue
        get(params, "tissue", nothing) != tissue && continue
        get(params, "modelid", nothing) != modelid && continue
        it = get(params, "it", nothing)
        isnothing(it) && continue
        files[it] = file
    end
    return files
end

## ---------------------------------------------------------------------
# stoi err
function _plot_stoi_err(modelid, tissue)

    files = _find_epout_files(modelid, tissue)

    model = AG.load_model(;modelid, tissue, uncompress = false)
    M, N = size(model)

    # get fluxs
    cid = (:AVS_CACHE, modelid, tissue)
    # delcache(AG, cid) # Reset
    its, avs = lcache(AG, cid) do
        _avs = Vector{Vector{Float64}}()
        _its = Vector{Int}()
        for (it, file) in sort!(collect(files)) # Test
            _epout = ldat(file)
            push!(_its, it)
            push!(_avs, MetEP.av(_epout))
        end
        return _its, _avs
    end

    # sub range
    itrange = 1:80
    its = its[itrange]
    avs = avs[itrange]

    cid = (:ERRS_CACHE, modelid, tissue, its)
    # delcache(AG, cid) # Reset
    abs_errs, norm_errs = lcache(AG, cid) do
        _abs_errs, _norm_errs = [], []
        @showprogress for (it, av) in zip(its, avs)
            _abs_err = MetNets.stoi_err(model, av)
            _norm_err = MetNets.norm_stoi_err(model, av; normfun = (metv) -> maximum(abs, metv))
            push!(_abs_errs, _abs_err)
            push!(_norm_errs, _norm_err)
        end
        return _abs_errs, _norm_errs
    end

    # get errs
    stoi_errav, stoi_errstd = Float64[], Float64[]
    stoi_norm_errav, stoi_norm_errstd = Float64[], Float64[]
    stoi_norm2_errav, stoi_norm2_errstd = Float64[], Float64[]
    errs_samples = Dict{Int, Vector{Float64}}()
    samples_its = its[1:end ÷ 4:end]
    @show samples_its
    @showprogress for (it, abs_err, norm_err) in zip(its, abs_errs, norm_errs)

        # trend
        
        push!(stoi_errav, maximum(abs, abs_err))
        push!(stoi_errstd, std(abs_err))

        push!(stoi_norm_errav, maximum(abs, norm_err))
        push!(stoi_norm_errstd, std(norm_err))
        
        push!(stoi_norm2_errav, mean(abs, norm_err))
        push!(stoi_norm2_errstd, 0.0)
        
        # hist
        if it in samples_its
            errs_samples[it] = abs.(norm_err)
        end

    end

    # Err trend
    _fun(arr) = log10.(arr)
    p1 = plot(its, _fun(stoi_errav); 
        label = "", lw = 3, c = :black, xlabel = "EP it", ylabel = "abs stoi err", title = modelid, 
        ylim = [0.0, Inf]
    )
    plot!(p1, its, _fun(stoi_errav .- stoi_errstd); label = "", lw = 1, ls = :dash, c = :black)
    plot!(p1, its, _fun(stoi_errav .+ stoi_errstd); label = "", lw = 1, ls = :dash, c = :black)
    
    p2 = plot(its, _fun(stoi_norm_errav); 
        label = "", lw = 3, c = :black, xlabel = "EP it", ylabel = "maximum norm stoi err", title = modelid,
        ylim = [0.0, Inf]
    )
    plot!(p2, its, _fun(stoi_norm_errav .- stoi_norm_errstd); label = "", lw = 1, ls = :dash, c = :black)
    plot!(p2, its, _fun(stoi_norm_errav .+ stoi_norm_errstd); label = "", lw = 1, ls = :dash, c = :black)

    p3 = plot(its, _fun(stoi_norm2_errav); 
        label = "", lw = 3, c = :black, xlabel = "EP it", ylabel = "mean norm stoi err", title = modelid
    )
    plot!(p3, its, _fun(stoi_norm2_errav .- stoi_norm2_errstd); label = "", lw = 1, ls = :dash, c = :black)
    plot!(p3, its, _fun(stoi_norm2_errav .+ stoi_norm2_errstd); label = "", lw = 1, ls = :dash, c = :black)
    
    ps = [p1, p2, p3]
    sfig(AG, ps,
        @fileid, "maxent_ep_beta0", "stoi_err", (;modelid, tissue), ".png";
        layout = (1, length(ps))
    )

    # Err hist
    ps = Plots.Plot[]
    for (it, errs) in sort!(collect(errs_samples))

        # Show big errs
        if it == last(samples_its)
            _max = maximum(abs, errs)
            _th = 0.1
            _big_erris = findall(errs .> _max * _th)
            
            println()
            println("-"^54)
            println("-"^54)
            @info("Doing", it, errth =  _max * _th)
            for meti in _big_erris
                println("-"^54)
                MetNets.summary(model, model.mets[meti])
                println()
                println("rel stoi err: ", errs[meti])
                println()
            end
            println()
        end
        
        p = histogram(_fun(errs); 
            label = "", title = string(modelid, " it: ", it), 
            bins = 200, xlabel = "stoi err"
        )
        push!(ps, p)
    end

    sfig(AG, ps,
        @fileid, "maxent_ep_beta0", "stoi_err_errs_samples", (;modelid, tissue), ".png";
        layout = (1, length(ps))
    )

    # Err vs met_rxns
    lit = last(samples_its)
    lerrs = errs_samples[lit]
    p = plot(; xlabel = "log rxn count", ylabel = "log stoi err", title = string(modelid, " it: ", lit))
    rxns_counts = Int[]
    _errs = Float64[]
    for (meti, err) in enumerate(lerrs)
        rxnis = MetNets.met_rxns(model, meti)
        rxns_count = length(rxnis)
        rxns_count == 0 && continue
        err == 0 && continue
        if rxns_count > 100
            MetNets.summary(model, model.mets[meti])
        end
        push!(rxns_counts, rxns_count)
        push!(_errs, err)
    end
    scatter!(p, log10.(rxns_counts), log10.(_errs); label = "", c = :black)
    sfig(AG, p,
        @fileid, "maxent_ep_beta0", "stoi_err_vs_rxns_count", (;modelid, tissue, it = lit), ".png"
    )

end

let
    tissue = "GBM"
    for modelid in ["fva_base", "fva_scaled"]
        _plot_stoi_err(modelid, tissue)
    end
end

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    modelid = "base"
    # model = AG.load_model(;modelid, tissue, uncompress = false) 
    model = HG.load_humangem_base_model(;uncompress = false)
    
    big_err_mets = [ "m01401s", "m02394c", "m02394s", "m02680c", "m02680s", "m02871m", "m02877m" ]
    for met in big_err_mets
        println("-"^54)
        MetNets.summary(model, met)
        println()
    end
    # MetNets.balance_str.([model], big_err_mets)
end

## ---------------------------------------------------------------------
# let
#     tissue = "GBM"
#     modelid = "fva_base"
#     model = AG.load_model(modelid, tissue; uncompress = false) 
#     MetNets.similar_rxns(model)
# end