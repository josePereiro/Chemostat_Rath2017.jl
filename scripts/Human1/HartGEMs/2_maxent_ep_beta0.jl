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
# DESCRIPTION
# This is a numeric test of MaxEnt EP at beta 0

## ---------------------------------------------------------------------
let
    return
    tissue = "GBM"
    for modelid in ["fva_base", "fva_scaled"]

        model = AG.load_model(modelid, tissue; uncompress = true)
        model = MetNets.force_dims(model)
        objider = HG.HUMAN_BIOMASS_IDER
        M, N = size(model)
        objidx = MetNets.rxnindex(model, objider)

        # maxent
        save_frec = 5
        function oniter(it, epmodel)
            !iszero(rem(it, save_frec)) && return
            
            epout_ = MetEP.produce_epout(epmodel; drop_epfields = true)
            sdat(AG, epout_,
                "maxent_ep_beta0", "epout", (;modelid, tissue, it), ".jls";
                verbose = false
            )
            
        end
        
        try
            @info("MAxEnt", model = size(model))
            epmodel = MetEP.EPModel(model)
            model = nothing; GC.gc(); # free mem
            MetEP.converge_ep!(epmodel; oniter,
                maxiter = 2000, # Test 
                verbose = true
            )
        catch err
            @error("At MaxEnt ", err)
        end

        # dat = (;stoi_errav, stoi_norm_errav, stoi_errstd, stoi_norm_errstd)
        # sdat(AG, dat, 
        #     "stoi_errs", (;modelid, tissue), ".jls";
        #     verbose = true
        # )
    end

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

    model = AG.load_model(modelid, tissue; uncompress = false)

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
    range = 1:60
    its = its[range]
    avs = avs[range]

    # get errs
    stoi_errav, stoi_errstd = Float64[], Float64[]
    stoi_norm_errav, stoi_norm_errstd = Float64[], Float64[]
    @showprogress for (it, av) in zip(its, avs)

        errs = MetNets.stoi_err(model, av)
        push!(stoi_errav, maximum(abs, errs))
        push!(stoi_errstd, std(errs))

        errs = MetNets.norm_stoi_err(model, av; normfun = (metv) -> mean(abs, metv))
        push!(stoi_norm_errav, maximum(abs, errs))
        push!(stoi_norm_errstd, std(errs))

    end

    # Plot
    _fun(arr) = arr
    p1 = plot(its, _fun(stoi_errav); label = "", lw = 3, c = :black, xlabel = "EP it", ylabel = "abs stoi err", title = modelid)
    plot!(p1, its, _fun(stoi_errav .- stoi_errstd); label = "", lw = 1, ls = :dash, c = :black)
    plot!(p1, its, _fun(stoi_errav .+ stoi_errstd); label = "", lw = 1, ls = :dash, c = :black)
    
    p2 = plot(its, _fun(stoi_norm_errav); label = "", lw = 3, c = :black, xlabel = "EP it", ylabel = "maximum norm stoi err", title = modelid)
    plot!(p2, its, _fun(stoi_norm_errav .- stoi_norm_errstd); label = "", lw = 1, ls = :dash, c = :black)
    plot!(p2, its, _fun(stoi_norm_errav .+ stoi_norm_errstd); label = "", lw = 1, ls = :dash, c = :black)
    
    sfig(AG, [p1, p2],
        @fileid, "maxent_ep_beta0_stoi_err", (;modelid, tissue), ".png"
    )

end

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    for modelid in ["fva_base", "fva_scaled"]
        _plot_stoi_err(modelid, tissue)
    end
end

## ---------------------------------------------------------------------
# ## ---------------------------------------------------------------------
# # MAXENT
# let
#     tissue = "GBM"

#     for modelid in ["fva_base", "fva_scaled"]
#         (modelid == "fva_base") && continue

#         # load model
#         model = AG.load_model(modelid, tissue; uncompress = true)
        
#         # MaxEnt
#         @info("Doing", tissue, modelid, model = size(model))
#         conv_err = Float64[]
#         function oniter(_, epmodel)
#             err = get(epmodel.stat, :max_err, nothing)
#             isnothing(err) && return
#             push!(conv_err, err)
#             return false, nothing
#         end

#         # EP
#         alpha = Inf
#         epsconv = 1e-4
#         try
#             epout = Ch.MaxEntEP.maxent_ep(model;
#                 alpha,
#                 verbose = true, 
#                 epsconv,
#                 maxiter = 1000, 
#                 oniter
#             )
#         catch err
#             sdat(AG, (;conv_err, alpha, epsconv), 
#                 "maxent_ep_beta0_conv_err", (;modelid, tissue), ".jls";
#                 verbose = true 
#             )
#             @error("ERROR", err)
#             continue
#         end
#         ep_z = Ch.Utils.av(model, epout, HG.HUMAN_BIOMASS_IDER)

#         # FBA
#         fbaout = MetLP.fba(model, HG.HUMAN_BIOMASS_IDER);
#         fba_z = Ch.Utils.av(model, fbaout, HG.HUMAN_BIOMASS_IDER)

#         # Exp
#         exp_z = maximum(Rd.val(:μ))

#         println("\n\n")
#         @info("Done", tissue, ep_z, fba_z, exp_z)
#     end
# end

# ## ---------------------------------------------------------------------
# using Plots

# ## ---------------------------------------------------------------------
# let
#     tissue = "GBM"
#     ps = Plots.Plot[]
#     for modelid in ["fva_base", "fva_scaled"]
#         base_model = AG.load_model(modelid, tissue; uncompress = false)
#         @show size(base_model)
#         nzS = filter((s) -> !iszero(s), base_model.S[:])
#         absS = abs.(nzS)
#         p = histogram(log10.(absS); 
#             label = "", title = "$(tissue) $(modelid)", 
#             xlabel = "log10(S)", yaxis = nothing
#         )
#         push!(ps, p)
#     end
    
#     sfig(AG, ps,
#         @fileid, "maxent_ep_beta0_S_hist", ".png"
#     )
# end;

# ## ---------------------------------------------------------------------
# let
#     tissue = "GBM"
#     ps = map(["fva_base", "fva_scaled"]) do modelid
#         conv_err, alpha, epsconv = ldat(AG, 
#             "maxent_ep_beta0_conv_err", (;modelid, tissue), ".jls"
#         )

#         conv_err = conv_err[1:end-1]
#         p = plot(log10.(conv_err); lw = 3, title = string(modelid),
#             label = "ep", xlabel = "iteration", ylabel = "log(converr)"
#         )

#         plot!(eachindex(conv_err), fill(log10(epsconv), length(conv_err));
#             label = "target", ls = :dash, lw = 3
#         )
#     end
#     sfig(AG, ps,
#         @fileid, "maxent_ep_beta0_conv_err", ".png"
#     )
# end

# ## ---------------------------------------------------------------------
# using LinearAlgebra
# let
#     tissue = "GBM"
#     ps = map(["fva_base", "fva_scaled"]) do modelid
#         model = AG.load_model(modelid, tissue; uncompress = false)
#         p = histogram(abs.(model.ub - model.lb); 
#             label = "", title = "$(tissue) $(modelid)", 
#             xlabel = "reaction", yaxis = "ub - lb"
#         )
#     end
#     sfig(AG, ps,
#         @fileid, "degeneracy", ".png"
#     )
# end