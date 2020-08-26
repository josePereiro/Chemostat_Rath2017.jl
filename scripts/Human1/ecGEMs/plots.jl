## Ploting

using Plots
pyplot();

ider_map = Dict(rath_ider => HG.exch_met_map[HG.mets_map[rath_ider]] for rath_ider in Rd.msd_mets)
ider_map["μ"] = OBJ_IDER
for (id1, id2) in ider_map
    ider_map[id2] = id1
end


## FBA vs ξ
p = plot(title = "FBA", xlabel = "ξ", ylabel = "flx", xscale = :log10)
for (model_id, model) in ec_models
    for (stst, fbaouts) in fbaouts_dict[model_id]
        growth_avs = []
        cost_avs = []
        for (ξ, fbaout) in zip(ξs, fbaouts)
            growth = Ch.Utils.av(model, fbaout, OBJ_IDER)
            cost = Ch.Utils.av(model, fbaout, PROT_POOL_EXCHANGE)

            push!(growth_avs, growth)
            push!(cost_avs, cost)
        end
        plot!(p, ξs, growth_avs, label = "")
        # plot!(p, ξs, cost_avs, label = "")
    end
end
scatter!(p, Rd.val(:ξ, Rd.ststs), Rd.val(:μ, Rd.ststs), label = "")
p

## Correlations
lim = 0.2
lims = [-lim, lim]
p = plot(
    ylim = lims,
    xlim = lims,
    title = "FBA correlation",
    xlabel = "exp",
    ylabel = "model"
)
for (model_id, ec_model) in ec_models
    for (stst, fbaouts) in fbaouts_dict[model_id]
        exp_ξ =  Rd.val(:ξ, stst)
        # exp_ξ
        exp_ξi = findfirst(isequal(exp_ξ), ξs)
        fbaout = fbaouts[exp_ξi]
        for rider in [Rd.msd_mets; "μ"]
            sense = rider == "μ" ? 1 : -1
            mider = ider_map[rider]
            exp_av = sense * Rd.qval(rider, stst)
            exp_err = Rd.qerr(rider, stst)
            mod_av = Ch.Utils.av(ec_model, fbaout, mider)
            scatter!(p, [exp_av], [mod_av]; xerr = [exp_err], label = "")
        end
    end
end
p

## Test
