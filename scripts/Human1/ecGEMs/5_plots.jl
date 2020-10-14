## ---------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## ---------------------------------------------------------------
import Plots: plot, plot!, scatter, scatter!
import SparseArrays
import Distributions: mean
import Chemostat_Rath2017: Chemostat
import Chemostat_Rath2017.Chemostat.LP: MathProgBase
import Chemostat_Rath2017.Human1: HumanGEM, RathData, ecGEMs
const ChU = Chemostat.Utils
const HG = HumanGEM
const Rd = RathData
const ecG = ecGEMs

## ---------------------------------------------------------------
maxent_dat = ChU.load_data(ecG.MAXENT_FBA_EB_BOUNDLES_FILE)

## ---------------------------------------------------------------
commit_hash = ChU.load_commit_short_hash(ecG.MAXENT_FBA_EB_BOUNDLES_FILE)
fig_title = string(nameof(ecG), " [", commit_hash, "]")
println(fig_title)

## ------------------------------------------------------------------
color_pool = [:orange, :blue, :red, :black, :violet, 
    :gray, :green, :brown, :magenta];
colors = Dict(exp => rand(color_pool) for exp in Rd.exps)

## ------------------------------------------------------------------
closest_βs = Dict()
for (model_id, bundles) in maxent_dat
    closest_βs[model_id] = Dict()
    for (exp, bundle) in bundles
        exp_xi = Rd.val("ξ", exp)
        exp_mu = Rd.val("μ", exp)
        closest_βs[model_id][exp] = ChU.find_closest_beta(bundle, exp_xi, exp_mu, HG.OBJ_IDER)
    end
end
closest_βs

## ---------------------------------------------------------------
# GROWTH vs BETA
function growth_vs_beta()
    eps_ = 1e-5
    p = plot(title = fig_title, 
        xlabel = "log beta", ylabel = "log growth rate", 
        legend = false
    )
    serie_fun = log10 #(x) -> x
    for (model_id, bundles) in maxent_dat
        sbundles = sort(collect(bundles); by = first)
        for (exp, bundle) in sbundles
            # exp != 7 && continue
            exp_growth = Rd.val("μ", exp)
            plot!(p, serie_fun.(bundle.βs), fill(serie_fun(exp_growth), length(bundle.βs)); 
                ls = :dash, color = colors[exp], lw = 3, label = "")
            
            exp_xi = Rd.val("ξ", exp)
            fba_growth = ChU.av(bundle, exp_xi, :fba, HG.OBJ_IDER)
            fba_growth += eps_
            plot!(p, serie_fun.(bundle.βs), fill(serie_fun(fba_growth), length(bundle.βs)); 
                ls = :dot, color = colors[exp], lw = 3, label = "")

            ep_growths = ChU.av(bundle, exp_xi, bundle.βs, :ep, HG.OBJ_IDER)
            ep_growths .+= eps_
            ep_stds = sqrt.(ChU.va(bundle, exp_xi, bundle.βs, :ep, HG.OBJ_IDER))
            plot!(p, serie_fun.(bundle.βs), serie_fun.(ep_growths); 
                label = exp, lw = 3, color = colors[exp])

        end
    end
    return p
end
growth_vs_beta()

## ------------------------------------------------------------------
# TOTAL STOI ERROR vs BETA
function total_stoi_err_vs_beta()
    p = plot(title = fig_title, 
        xlabel = "log beta", ylabel = "stoi error [min/mean/max]")
        
    for (model, bundles) in maxent_dat

        sbundles = sort(collect(bundles); by = (x) -> x[1])
        for (exp, bundle) in sbundles
            
            exp_xi = Rd.val("ξ", exp)
            model = bundle[exp_xi, :net]

            # collect errors
            M, N = size(model)
            errs = zeros(M, length(bundle.βs))
            for (βi, β) in bundle.βs |> enumerate
                epout = bundle[exp_xi, β, :ep]
                errs[:, βi] = ChU.norm1_stoi_err(model, epout)
            end

            # plots
            plot!(p, log10.(bundle.βs), maximum.(eachcol(errs));
                lw = 3, ls = :dash, c = colors[exp], label = "")
            plot!(p, log10.(bundle.βs), mean.(eachcol(errs));
                lw = 3, ls = :solid , c = colors[exp], label = string(exp))
            plot!(p, log10.(bundle.βs), minimum.(eachcol(errs));
                lw = 3, ls = :dot, c = colors[exp], label = "")
        end
    end
    return p
end
total_stoi_err_vs_beta()

## ------------------------------------------------------------------
# METS STOI ERROR vs BETA
function mets_stoi_err_vs_beta(model_id = "GTEx-brain", exp = "E")

    exch_met_map = HG.load_exch_met_map()
    p = plot(title = string(fig_title, " ", exp), 
            xlabel = "log beta", ylabel = "stoi error ")

    bundle = maxent_dat[model_id][exp]
    exp_xi = Rd.val("ξ", exp)
    model = bundle[exp_xi, :net]
    idermap = HG.load_readable_met_ids_map();

    # collect errors
    err_dict = Dict()
    for (exch, _) in model.intake_info
        met = exch_met_map[exch]
        for (βi, β) in bundle.βs |> enumerate
            epout = bundle[exp_xi, β, :ep]
            errs = get!(err_dict, met, [])
            push!(errs, ChU.norm1_stoi_err(model, epout, met))
        end
    end

    # plots
    serr_dict = sort(collect(err_dict); by = (p) -> sum(p[2]), rev = true)
    lcount = 5
    serie_fun = log10
    for (ider, errs) in serr_dict
        lb_ = lcount < 0 ? "" : get(idermap, ider, ider)
        plot!(p, serie_fun.(bundle.βs), errs;
            lw = 3, ls = :dash, label = lb_)
        lcount -= 1
    end
    return p
end
mets_stoi_err_vs_beta()

## -------------------------------------------------------------------
# CLOSEST BETA VS BETA
# Just for checking that the experimental objective is inside the beta intervals
# and evaluate the 'experimental' beta approximation
function closest_beta()
    p = plot(title = fig_title, 
            xlabel = "log beta", ylabel = "growth rate"
    )
    
    serie_fun = log10
    for (model_id, bundles) in maxent_dat

        sbundles = sort(collect(bundles); by = (x) -> x[1])
        for (exp, bundle) in sbundles
            
            # model
            exp_xi = Rd.val("ξ", exp)
            μs = ChU.av(bundle, exp_xi, bundle.βs, :ep, HG.OBJ_IDER)
            plot!(p, serie_fun.(bundle.βs), μs, label = "", color = colors[exp], lw = 3)

            # exp
            scatter!(p, [serie_fun(closest_βs[model_id][exp])], [Rd.val("μ", exp)], label = "", 
                color = colors[exp], ms = 9)

        end
    end
    return p
end
closest_beta();

## -------------------------------------------------------------------
# GROWTH CORRELATION
function growth_correlation()
    p = plot(title = fig_title, 
        xlabel = "experimental growth rate", ylabel = "modeled growth rate"
    )

    ider_fun(x) = x#log10(x)
    for (model_id, bundles) in maxent_dat

        sbundles = sort(collect(bundles); by = (x) -> x[1])
        for (exp, bundle) in sbundles
        
            exp_ξ =  Rd.val("ξ", exp)
            exp_β = closest_βs[model_id][exp]
            model_μ = ChU.av(bundle, exp_ξ, exp_β, :ep, HG.OBJ_IDER)
            exp_μ = Rd.val("μ", exp)
            scatter!(p, [ider_fun(exp_μ)], [ider_fun(model_μ)], label = "", color = colors[exp])
        end
    end
    return p
end
growth_correlation()

## ------------------------------------------------------------------
# TOTAL CORRELATIONS
function total_correlation()
    ider_map = HG.load_rath_met_exch_map()
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
    ticks = round.(collect(range(-1, 5; length = 4)), digits = 2)
    lims = [minimum(ticks), maximum(ticks)]
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
    for (model_id, bundles) in maxent_dat

        sbundles = sort(collect(bundles); by = (x) -> x[1])
        for (stst, bundle) in sbundles

            exp_ξ =  Rd.val(:ξ, stst)
            exp_μ =  Rd.val(:D, stst)
            exp_β = closest_βs[model_id][exp]

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
                    println("FBA\\", model_id, "\\", stst, "\\", rider, 
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
                     println("EP\\", model_id, "\\", stst, "\\", rider, 
                        ":\t [", rf(exp_av), ", ", rf(mod_av), ", ", rf(rel_err), "]")
                 end
                 println()

            end
        end
    end

    # FBA
    scatter!(fbap, xs, fbays; xerr = xerrs, label = "", color = :black)

    # EP
    scatter!(epp, xs, epys; xerr = xerrs, yerr = epyerrs, label = "", color = :black)

    return plot([fbap, epp]..., layout = 2, size = [800, 400],)
end
total_correlation()

## ------------------------------------------------------------------
# MARGINALS
function exp_b_marginals(;model_id = "GTEx-brain", exp = "E",
        xlims = Dict())
    bundle = maxent_dat[model_id][exp]
    exp_ξ =  Rd.val(:ξ, exp)
    exp_β = closest_βs[model_id][exp]
    model = bundle[exp_ξ, :net]
    ider_map = HG.load_rath_met_exch_map()

    ps = []
    for rider in Rd.iders_to_plot
        mider = ider_map[rider]
        p = plot(title = first(rider, 10), 
            titlefont = 8, 
            yticks = []
        )
        ChP.plot_marginal!(p, bundle, exp_ξ, exp_β, [:fba, :ep], mider;
             label = "", xlim = get(xlims, rider, nothing))
        push!(ps, p)
    end
    plot(ps..., layout = length(ps))
end
exp_b_marginals(exp = "E"; 
    xlims = Dict("GLC" => [-5, 5], 
                "LAC" => [-10, 250],
                "GLN" => [-1, 1],
                "NH4" => [-1, 8],
                "GAL" => [-1, 10],
                "PYR" => [-10, 600],
                "GLU" => [-5, 5],
                "ALA" => [-5, 50],
                "ASP" => [-5, 5],
                "μ" => [-0.1, 0.1], 
    )
)
