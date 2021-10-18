using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using ProgressMeter
    using Base.Threads

    import Chemostat
    const Ch = Chemostat
    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
    const AG = H1.HartGEMs
end

## ---------------------------------------------------------------------
# DESCRIPTION
# This is a numeric test of MaxEnt EP at beta 0

## ---------------------------------------------------------------------
# MAXENT
let
    stst = "A"
    tissue = "GBM"

    for modelid in ["fva_base", "fva_scaled"]
        (modelid == "fva_base") && continue

        # load model
        model = AG.load_model(modelid, tissue; uncompress = true)

        # Set medium
        ξ = Rd.val(:ξ, stst)
        base_intake_info = HG.stst_base_intake_info(stst) 
        Chemostat.apply_bound!(model, ξ, base_intake_info; 
            emptyfirst = true, ignore_miss = true    
        )
        
        # MaxEnt
        @info("Doing", stst, tissue, modelid, model = size(model))
        conv_err = Float64[]
        function oniter(_, epmodel)
            err = get(epmodel.stat, :max_err, nothing)
            isnothing(err) && return
            push!(conv_err, err)
            return false, nothing
        end

        # EP
        alpha = Inf
        epsconv = 1e-4
        try
            epout = Ch.MaxEntEP.maxent_ep(model;
                alpha,
                verbose = true, 
                epsconv,
                maxiter = 1000, 
                oniter
            )
        catch err
            sdat(AG, (;conv_err, alpha, epsconv), 
                "maxent_ep_beta0_conv_err", (;modelid, tissue), ".jls";
                verbose = true 
            )
            @error("ERROR", err)
            continue
        end
        ep_z = Ch.Utils.av(model, epout, HG.HUMAN_BIOMASS_IDER)

        # FBA
        fbaout = Chemostat.LP.fba(model, HG.HUMAN_BIOMASS_IDER);
        fba_z = Ch.Utils.av(model, fbaout, HG.HUMAN_BIOMASS_IDER)

        # Exp
        exp_z = Rd.val(:μ, stst)

        println("\n\n")
        @info("Done", stst, tissue, ep_z, fba_z, exp_z)
    end
end

## ---------------------------------------------------------------------
using Plots

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    ps = Plots.Plot[]
    for modelid in ["fva_base", "fva_scaled"]
        base_model = AG.load_model(modelid, tissue; uncompress = false)
        @show size(base_model)
        nzS = filter((s) -> !iszero(s), base_model.S[:])
        absS = abs.(nzS)
        p = histogram(log10.(absS); 
            label = "", title = "$(tissue) $(modelid)", 
            xlabel = "log10(S)", yaxis = nothing
        )
        push!(ps, p)
    end
    
    sfig(AG, ps,
        @fileid, "maxent_ep_beta0_S_hist", ".png"
    )
end;

## ---------------------------------------------------------------------
let
    tissue = "GBM"
    ps = map(["fva_base", "fva_scaled"]) do modelid
        conv_err, alpha, epsconv = ldat(AG, 
            "maxent_ep_beta0_conv_err", (;modelid, tissue), ".jls"
        )

        conv_err = conv_err[1:end-1]
        p = plot(log10.(conv_err); lw = 3, title = string(modelid),
            label = "ep", xlabel = "iteration", ylabel = "log(converr)"
        )

        plot!(eachindex(conv_err), fill(log10(epsconv), length(conv_err));
            label = "target", ls = :dash, lw = 3
        )
    end
    sfig(AG, ps,
        @fileid, "maxent_ep_beta0_conv_err", ".png"
    )
end

## ---------------------------------------------------------------------
using LinearAlgebra
let
    tissue = "GBM"
    ps = map(["fva_base", "fva_scaled"]) do modelid
        model = AG.load_model(modelid, tissue; uncompress = false)
        p = histogram(abs.(model.ub - model.lb); 
            label = "", title = "$(tissue) $(modelid)", 
            xlabel = "reaction", yaxis = "ub - lb"
        )
    end
    sfig(AG, ps,
        @fileid, "degeneracy", ".png"
    )
end