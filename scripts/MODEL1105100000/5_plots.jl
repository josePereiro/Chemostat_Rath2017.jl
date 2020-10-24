import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

import Plots: plot, plot!, scatter, scatter!
import SparseArrays
import Distributions: mean
import Chemostat_Rath2017: Chemostat, MODEL1105100000
import Chemostat_Rath2017.Chemostat.LP: MathProgBase
import Chemostat_Rath2017.Human1: HumanGEM, RathData
const ChU = Chemostat.Utils
const ChP = Chemostat.Plots
const HG = HumanGEM
const Rd = RathData
const M = MODEL1105100000

## ------------------------------------------------------------------
# LOADING BUNDLES
bundles = ChU.load_data(M.MAXENT_FBA_EB_BUNDLES_FILE)

## ---------------------------------------------------------------
commit_hash = ChU.load_commit_short_hash(M.MAXENT_FBA_EB_BUNDLES_FILE)
fig_title = string(nameof(M), " [", commit_hash, "]")
println(fig_title)

## ------------------------------------------------------------------
color_pool = [:orange, :blue, :red, :black, :violet, 
    :gray, :green, :brown, :magenta];
colors = Dict(exp => rand(color_pool) for exp in Rd.exps)

## ------------------------------------------------------------------
closest_βs = Dict()
for (exp, bundle) in bundles
    exp_xi = Rd.val("ξ", exp)
    exp_mu = Rd.val("μ", exp)
    closest_βs[exp] = ChU.find_closest_beta(bundle, exp_xi, exp_mu, M.BIOMASS_IDER)
end
closest_βs

## ------------------------------------------------------------------
function get_limits(xdat, ydat; marginf = 0.10)
    xmargin = abs(minimum(xdat) - maximum(xdat)) * marginf
    xlims = [ minimum(xdat) - xmargin, maximum(xdat) + xmargin ]

    ymargin = abs(minimum(ydat) - maximum(ydat)) * marginf
    ylims = [ minimum(ydat) - ymargin, maximum(ydat) + ymargin ]

    return (xlims, ylims)
end

## ------------------------------------------------------------------
# TOTAL CORRELATION
function total_correlation()
    ider_map = M.load_rath_met_exch_map()
    ticks = round.(collect(range(-0.4, 0.4; length = 4)), digits = 2)
    lims = [minimum(ticks), maximum(ticks)]
    fbap = plot(
        ylim = lims, 
        xlim = lims, 
        xticks = ticks,
        yticks = ticks,
        xtickfontsize = 9,
        ytickfontsize = 9,
        title = "FBA correlation", 
        xlabel = "exp", ylabel = "model"
    )
    plot!(fbap, x -> x, lims[1], lims[2]; label = "", color = :black)
    # ticks = round.(collect(range(-0.4, 0.4; length = 4)), digits = 2)
    # lims = [minimum(ticks), maximum(ticks)]
    epp = plot(
        ylim = lims, 
        xlim = lims, 
        xticks = ticks,
        yticks = ticks,
        xtickfontsize = 9,
        ytickfontsize = 9,
        title = "EP correlation", 
        xlabel = "exp", ylabel = "model"
    )
    plot!(epp, x -> x, lims[1], lims[2]; label = "", color = :black)

    # Info print
    rel_th = 0.2;
    abs_th = 1e-4;
    rf(x) = round(x, digits = 4)

    xs = []
    xerrs = []
    fbays = []
    epys = []
    epyerrs = []
    for (stst, bundle) in bundles
        exp_ξ =  Rd.val(:ξ, stst)
        exp_μ =  Rd.val(:D, stst)
        exp_β = closest_βs[stst]

        for rider in Rd.iders_to_plot
            sense = rider == Rd.growth_ider ? 1 : -1 # TODO: package this
            mider = ider_map[rider]
            exp_av = sense * Rd.qval(rider, stst)
            exp_err = Rd.qerr(rider, stst)
            push!(xs, exp_av)
            push!(xerrs, exp_err)

            # FBA
            mod_av = ChU.av(bundle, exp_ξ, :fba, mider)
            push!(fbays, mod_av)

            # Info print
            # exp\\mod\\rel
            rel_err = abs(mod_av - exp_av)/max(abs(mod_av), abs(exp_av))
            if abs(mod_av) > abs_th && rel_th < rel_err < Inf
                println("FBA", "\\", stst, "\\", rider, 
                    ":\t [", rf(exp_av), ", ", rf(mod_av), ", ", rf(rel_err), "]")
            end
            
            # EP
            mod_av = ChU.av(bundle, exp_ξ, exp_β,:ep, mider)
            push!(epys, mod_av)
            mod_err = sqrt(ChU.va(bundle, exp_ξ, exp_β,:ep, mider))
            push!(epyerrs, mod_err)
            
            # Info print
            rel_err = abs(mod_av - exp_av)/max(abs(mod_av), abs(exp_av))
            if abs(mod_av) > abs_th && rel_th < rel_err < Inf
                println("EP", "\\", stst, "\\", rider, 
                    ":\t [", rf(exp_av), ", ", rf(mod_av), ", ", rf(rel_err), "]")
            end
        end
    end

    # FBA
    scatter!(fbap, xs, fbays; xerr = xerrs, label = "", color = :black)

    # EP
    scatter!(epp, xs, epys; xerr = xerrs, yerr = epyerrs, label = "", color = :black)

    return plot([fbap, epp]...; title = fig_title, layout = 2, size = [800, 400],)
end
total_correlation()

## ------------------------------------------------------------------
# MARGINALS
function exp_b_marginals(exp = "E")
    bundle = bundles[exp]
    exp_ξ =  Rd.val(:ξ, exp)
    exp_β = closest_βs[exp]
    model = bundle[exp_ξ, :net]
    ider_map = M.load_rath_met_exch_map()

    ps = []
    for rider in Rd.iders_to_plot
        mider = ider_map[rider]
        p = plot(title = first(mider, 10), 
            titlefont = 8, 
            xticks = [], yticks = []
        )
        ChP.plot_marginal!(p, bundle, exp_ξ, exp_β, [:fba, :ep], mider, label = "")
        push!(ps, p)
    end
    plot(ps..., layout = length(ps))
end
exp_b_marginals("E")

## ------------------------------------------------------------------
# # INDIVIDUAL CORRELATIONS
# function individual_correlation()
    
#     fbaps =[]
#     epps = []
#     fontsize = 8
#     ntick = 4
#     size_ = [1000, 800]
#     for rider in Rd.iders_to_plot
#         fbap = plot(
#             xtickfontsize = fontsize,
#             ytickfontsize = fontsize,
#             titlefont = fontsize,
#             guidefontsize = fontsize,
#             title = "$rider FBA correlation", 
#             xlabel = "exp", ylabel = "model"
#         )

#         epp = plot(
#             xtickfontsize = fontsize,
#             ytickfontsize = fontsize,
#             titlefont = fontsize,
#             guidefontsize = fontsize,
#             title = "$rider EP correlation", 
#             xlabel = "exp", ylabel = "model"
#         )

#         exp_avs = []
#         exp_errs = []
#         fba_avs = []
#         ep_avs = []
#         ep_errs = []
#         for (stst, bundle) in bundles
#             exp_ξ =  Rd.val(:ξ, stst)
#             exp_μ =  Rd.val(:D, stst)
#             exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, M.BIOMASS_IDER)

#             sense = rider == Rd.growth_ider ? 1 : -1 # TODO: package this
#             mider = ider_map[rider]
#             exp_av = sense * Rd.qval(rider, stst)
#             exp_err = Rd.qerr(rider, stst)
#             push!(exp_avs, exp_av)
#             push!(exp_errs, exp_err)

#             # FBA
#             fba_av = av(bundle, exp_ξ, :fba, mider)
#             push!(fba_avs, fba_av)

#             # EP
#             ep_av = av(bundle, exp_ξ, exp_β,:ep, mider)
#             push!(ep_avs, ep_av)
#             ep_err = sqrt(va(bundle, exp_ξ, exp_β,:ep, mider))
#             push!(ep_errs, ep_err)
            
#         end

#         xlims, ylims = get_limits(exp_avs, fba_avs; marginf = 0.20)
#         scatter!(fbap, exp_avs, fba_avs; 
#             xlim = xlims, ylim = ylims,
#             xtick = round.(range(minimum(xlims), maximum(xlims); 
#                 length = ntick); digits = 3),
#             ytick =  round.(range(minimum(ylims), maximum(ylims); 
#                 length = ntick); digits = 3),
#             xerr = exp_errs, 
#             size = size_,
#             label = ""
#         )
#         plot!(fbap, x -> x, xlims[1], xlims[2]; label = "", color = :black)

#         xlims, ylims = get_limits(exp_avs, ep_avs; marginf = 0.20)
#         scatter!(epp, exp_avs, ep_avs; 
#             xlim = xlims, ylim = ylims,
#             xtick = round.(range(minimum(xlims), maximum(xlims); 
#                 length = ntick); digits = 3),
#             ytick =  round.(range(minimum(ylims), maximum(ylims); 
#                 length = ntick); digits = 3),
#             yerr = ep_errs, 
#             xerr = exp_errs, 
#             size = size_,
#             label = ""
#         )
#         plot!(epp, x -> x, xlims[1], xlims[2]; label = "", color = :black)
        
#         push!(fbaps, fbap)
#         push!(epps, epp)

#     end

#     # return plot([fbap, epp]..., size = [800, 400])
#     return fbaps, epps
# end
# fbaps, epps = individual_correlation();
# plot(epps...; layout = length(epps) )
# ##
# n = length(fbaps)
# plot(fbaps...; size = [n * 300, 300], layout = grid(1, n))

# ##



## ------------------------------------------------------------------
# Dev data
stst = "E"
bundle = bundles[stst]
exp_ξ = Rd.val(:ξ, stst)
exp_μ = Rd.val(:D, stst)
model = bundle[exp_ξ, :net]
fbaout = bundle[exp_ξ, :fba]
exp_β = ChU.find_closest_beta(bundle, exp_ξ, exp_μ, M.BIOMASS_IDER)
epout = bundle[exp_ξ, exp_β, :ep];