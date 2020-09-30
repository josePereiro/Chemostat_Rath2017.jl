using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

using SparseArrays
using MathProgBase
using Distributions

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.Utils: av, va, μ, σ, bounds, is_exchange,
                        norm1_stoi_err, norm2_stoi_err,
                        find_closest_beta, metindex, met_rxns,
                        load_data, to_symbol_dict
                        
import Chemostat_Rath2017
import Chemostat_Rath2017: RathData, Human1
import Chemostat_Rath2017.Human1: HumanGEM, ecGEMs
const H1 = Human1
const Rd = RathData
const HG = HumanGEM
const ecG = ecGEMs

## Loading data
src_file = ecG.MAXENT_FBA_EB_BOUNDLES_FILE
bundles = load_data(src_file)

## Dev data
model_id = "GTEx-brain"
stst = "E"
bundle = bundles[model_id][stst]
exp_ξ = Rd.val(:ξ, stst)
exp_μ = Rd.val(:D, stst)
model = bundle[exp_ξ, :net]
fbaout = bundle[exp_ξ, :fba]
exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, H1.OBJ_IDER)
epout = bundle[exp_ξ, exp_β, :ep];

##
using Plots
# gr(size = (600, 500))

## ------------------------------------------------------------------
for (model_id, dat) in bundles
    for (stst, bundle) in dat
        exp_ξ = Rd.val(:ξ, stst)
        exp_μ = Rd.val(:D, stst)
        exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, H1.OBJ_IDER)
        mod_μ = av(bundle, exp_ξ, exp_β, :ep, H1.OBJ_IDER)
        println(stst, "  \texp mu: ", exp_μ, "     \t mod mu: ", mod_μ, "        \tbeta: ", exp_β)
    end
end
## ------------------------------------------------------------------
function get_limits(xdat, ydat; marginf = 0.10)
    xmargin = abs(minimum(xdat) - maximum(xdat)) * marginf
    xlims = [ minimum(xdat) - xmargin, maximum(xdat) + xmargin ]

    ymargin = abs(minimum(ydat) - maximum(ydat)) * marginf
    ylims = [ minimum(ydat) - ymargin, maximum(ydat) + ymargin ]

    return xlims, ylims
end

## ------------------------------------------------------------------
# TOOLS
# Total Correlations
ider_map = HG.load_rath_met_exch_map()
function total_correlation()

    xs = []
    xerrs = []
    fbays = []
    epys = []
    epyerrs = []
    for (model_id, dat) in bundles
        for (stst, bundle) in dat
            exp_ξ =  Rd.val(:ξ, stst)
            exp_μ =  Rd.val(:D, stst)
            exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, HG.OBJ_IDER)

            for rider in Rd.iders_to_plot
                sense = rider == Rd.growth_ider ? 1 : -1 # TODO: package this
                mider = ider_map[rider]
                exp_av = sense * Rd.qval(rider, stst)
                exp_err = Rd.qerr(rider, stst)
                push!(xs, exp_av)
                push!(xerrs, exp_err)

                # FBA
                mod_av = av(bundle, exp_ξ, :fba, mider)
                push!(fbays, mod_av)
                
                # EP
                mod_av = av(bundle, exp_ξ, exp_β,:ep, mider)
                push!(epys, mod_av)
                mod_err = sqrt(va(bundle, exp_ξ, exp_β,:ep, mider))
                push!(epyerrs, mod_err)
            end
        end
    end

    ntick = 4
    ndig = 3
    font = 9
    
    # FBA
    xlims, ylims = get_limits(xs, fbays; marginf = 0.20)
    xticks = round.(collect(range(minimum(xlims), maximum(xlims); 
        length = ntick)), digits = ndig)
    yticks = round.(collect(range(minimum(ylims), maximum(ylims); 
        length = ntick)), digits = ndig)
    fbap = plot(title = "FBA correlation", 
        xlabel = "exp", ylabel = "model", titlefont = font,
        xtickfontsize = font, ytickfontsize = font)
    plot!(fbap, x -> x, xlims[1], xlims[2]; label = "", color = :black)
    zoom_fbap = deepcopy(fbap)
    scatter!(fbap, xs, fbays; xerr = xerrs, label = "", color = :black)
    scatter!(zoom_fbap, xs, fbays; xerr = xerrs, 
        xlims = [-0.05, 0.04],
        ylims = [-0.1, 0.1],
        label = "", color = :black)

    # EP
    xlims, ylims = get_limits(xs, epys; marginf = 0.20)
    yticks = round.(collect(range(minimum(ylims), maximum(ylims); 
        length = ntick)), digits = ndig)
    epp = plot(title = "EP correlation", 
        xlabel = "exp", ylabel = "model", titlefont = font,
        xtickfontsize = font, ytickfontsize = font)
    plot!(epp, x -> x, xlims[1], xlims[2]; label = "", color = :black)
    zoom_epp = deepcopy(epp)
    scatter!(epp, xs, epys; xerr = xerrs, yerr = epyerrs, 
        label = "", color = :black)
    scatter!(zoom_epp, xs, epys; xerr = xerrs, yerr = epyerrs, 
        xlims = [-0.05, 0.04],
        ylims = [-1, 12],
        label = "", color = :black)

    return plot([fbap, epp, zoom_fbap, zoom_epp]..., layout = 4)
end
total_correlation()

## ------------------------------------------------------------------
# individual correlations
function individual_correlation()
    
    fbaps =[]
    epps = []
    fontsize = 8
    ntick = 4
    size_ = [1000, 800]
    for rider in Rd.iders_to_plot
        fbap = plot(
            xtickfontsize = fontsize,
            ytickfontsize = fontsize,
            titlefont = fontsize,
            guidefontsize = fontsize,
            title = "$rider FBA correlation", 
            xlabel = "exp", ylabel = "model"
        )

        epp = plot(
            xtickfontsize = fontsize,
            ytickfontsize = fontsize,
            titlefont = fontsize,
            guidefontsize = fontsize,
            title = "$rider EP correlation", 
            xlabel = "exp", ylabel = "model"
        )

        exp_avs = []
        exp_errs = []
        fba_avs = []
        ep_avs = []
        ep_errs = []
        for (model_id, dat) in bundles
            for (stst, bundle) in dat
                exp_ξ =  Rd.val(:ξ, stst)
                exp_μ =  Rd.val(:D, stst)
                exp_β = find_closest_beta(bundle, exp_ξ, exp_μ, M.OBJ_IDER)

                sense = rider == Rd.growth_ider ? 1 : -1 # TODO: package this
                mider = ider_map[rider]
                exp_av = sense * Rd.qval(rider, stst)
                exp_err = Rd.qerr(rider, stst)
                push!(exp_avs, exp_av)
                push!(exp_errs, exp_err)

                # FBA
                fba_av = av(bundle, exp_ξ, :fba, mider)
                push!(fba_avs, fba_av)

                # EP
                ep_av = av(bundle, exp_ξ, exp_β,:ep, mider)
                push!(ep_avs, ep_av)
                ep_err = sqrt(va(bundle, exp_ξ, exp_β,:ep, mider))
                push!(ep_errs, ep_err)
                
            end
        end

        xlims, ylims = get_limits(exp_avs, fba_avs; marginf = 0.20)
        scatter!(fbap, exp_avs, fba_avs; 
            xlim = xlims, ylim = ylims,
            xtick = round.(range(minimum(xlims), maximum(xlims); 
                length = ntick); digits = 3),
            ytick =  round.(range(minimum(ylims), maximum(ylims); 
                length = ntick); digits = 3),
            xerr = exp_errs, 
            size = size_,
            label = ""
        )
        plot!(fbap, x -> x, xlims[1], xlims[2]; label = "", color = :black)

        xlims, ylims = get_limits(exp_avs, ep_avs; marginf = 0.20)
        scatter!(epp, exp_avs, ep_avs; 
            xlim = xlims, ylim = ylims,
            xtick = round.(range(minimum(xlims), maximum(xlims); 
                length = ntick); digits = 3),
            ytick =  round.(range(minimum(ylims), maximum(ylims); 
                length = ntick); digits = 3),
            yerr = ep_errs, 
            xerr = exp_errs, 
            size = size_,
            label = ""
        )
        plot!(epp, x -> x, xlims[1], xlims[2]; label = "", color = :black)
        
        push!(fbaps, fbap)
        push!(epps, epp)

    end

    # return plot([fbap, epp]..., size = [800, 400])
    return fbaps, epps
end
fbaps, epps = individual_correlation();
plot(epps...; layout = length(epps) )

## ------------------------------------------------------------------
# stoi_err_by_compartments
# errs = log.(norm2_stoi_err(model, epout) .+ 1e-3);
# histogram(errs; normalize = true, xlabel = "log err", ylabel = "prob dens")

function stoi_err_by_compartments(get_comp::Function, model, out)
    
    errs = norm2_stoi_err(model, out)

    ## Split error by compartment
    cerrs = Dict()
    for (meti, met) in model.mets |> enumerate
        c = get_comp(model, met) 
        push!(get!(cerrs, c, []), errs[meti])
    end

    ps = []
    xlims_ = [-6,4]
    bins_ = 100
    for (c, errs) in cerrs
        log_err = log10.(errs .+ 1e-5) 
        p = histogram(log_err; label = string(c), 
            xlabel = "log err", ylabel = "prob", 
            xlims = xlims_, bins = bins_,
            normalize = true)
        push!(ps, p)
    end
    ##
    plot(ps..., size = [1000, 400])

end

cerrs = stoi_err_by_compartments(model, epout) do model, met
    rxns = met_rxns(model, met)
    return any(is_exchange(model, rxn) for rxn in rxns) ? "ext" : "int"
end

##
errs = norm1_stoi_err(model, epout);
maximum

##
plot(["A", "B"], [2,3])