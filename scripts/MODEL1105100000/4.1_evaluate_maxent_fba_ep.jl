import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## ------------------------------------------------------------------
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

## ------------------------------------------------------------------
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
    closest_βs[exp] = ChU.find_closest_beta(bundle, exp_xi, exp_mu, M.OBJ_IDER)
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
    sbundles = sort(collect(bundles); by = first)
    for (exp, bundle) in sbundles
        # exp != "A" && continue
        exp_growth = Rd.val("μ", exp)
        plot!(p, serie_fun.(bundle.βs), fill(serie_fun(exp_growth), length(bundle.βs)); 
            ls = :dash, color = colors[exp], lw = 3, label = "")
        
        exp_xi = Rd.val("ξ", exp)
        fba_growth = ChU.av(bundle, exp_xi, :fba, M.OBJ_IDER)
        fba_growth += eps_
        plot!(p, serie_fun.(bundle.βs), fill(serie_fun(fba_growth), length(bundle.βs)); 
            ls = :dot, color = colors[exp], lw = 3, label = "")

        ep_growths = ChU.av(bundle, exp_xi, bundle.βs, :ep, M.OBJ_IDER)
        ep_growths .+= eps_
        ep_stds = sqrt.(ChU.va(bundle, exp_xi, bundle.βs, :ep, M.OBJ_IDER))
        plot!(p, serie_fun.(bundle.βs), serie_fun.(ep_growths); 
            label = exp, lw = 3, color = colors[exp])

    end
    return p
end
growth_vs_beta()

## ------------------------------------------------------------------
# TOTAL STOI ERROR vs BETA
function total_stoi_err_vs_beta()
    p = plot(title = fig_title, 
        xlabel = "log beta", ylabel = "stoi error [min/mean/max]")
        
    sbundles = sort(collect(bundles); by = first)
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
    return p
end
total_stoi_err_vs_beta()


## ------------------------------------------------------------------
# METS STOI ERROR vs BETA
function mets_stoi_err_vs_beta(exp = "E")

    exch_met_map = M.load_exch_met_map()
    p = plot(title = string(fig_title, " ", exp), 
            xlabel = "log beta", ylabel = "stoi error ")

    bundle = bundles[exp]
    exp_xi = Rd.val("ξ", exp)
    model = bundle[exp_xi, :net]
    idermap = M.load_rath_met_exch_map();

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
    serr_dict = sort(collect(err_dict); by = (p) -> sum(last(p)), rev = true)
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

    sbundles = sort(collect(bundles); by = first)
    for (exp, bundle) in sbundles
        
        # model
        exp_xi = Rd.val("ξ", exp)
        μs = ChU.av(bundle, exp_xi, bundle.βs, :ep, M.OBJ_IDER)
        plot!(p, serie_fun.(bundle.βs), μs, label = "", color = colors[exp], lw = 3)

        # exp
        scatter!(p, [serie_fun(closest_βs[exp])], [Rd.val("μ", exp)], label = "", 
            color = colors[exp], ms = 9)

    end
    return p
end
closest_beta()

## -------------------------------------------------------------------
# GROWTH CORRELATION
function growth_correlation()
    p = plot(title = fig_title, 
        xlabel = "experimental growth rate", ylabel = "modeled growth rate"
    )

    serie_fun(x) = x#log10(x)
    sbundles = sort(collect(bundles); by = first)
    for (exp, bundle) in sbundles
    
        exp_ξ =  Rd.val("ξ", exp)
        exp_β = closest_βs[exp]
        model_μ = ChU.av(bundle, exp_ξ, exp_β, :ep, M.OBJ_IDER)
        exp_μ = Rd.val("μ", exp)
        scatter!(p, [serie_fun(exp_μ)], [serie_fun(model_μ)], label = "", color = colors[exp])
    end
    return p
end
growth_correlation()

## ------------------------------------------------------------------
function to_csv(exp = "E"; sep = ",")
    csv = []
    line = ["rxn idx", "rxn id", "fba av", "ep av", "ep va", "lb", "ub", "rxn eq", "rxn subSys"]
    push!(csv, line)

    exp_ξ =  Rd.val("ξ", exp)
    exp_β = closest_βs[exp]
    bundle = bundles[exp]
    model = bundle[exp_ξ, :net]
    idermap = M.load_rath_met_exch_map();

    iders = collect(keys(model.intake_info))
    iders = [iders; [idermap[rider] for rider in Rd.iders_to_plot]]
    for rxn in unique(iders)
        line = []
        # rxn idx
        rxni = ChU.rxnindex(model, rxn)
        push!(line, rxni)
        # rxn id
        push!(line, rxn)
        # fba av
        push!(line, ChU.av(bundle, exp_ξ, :fba, rxni))
        # ep av
        push!(line, ChU.av(bundle, exp_ξ, exp_β, :ep, rxni))
        # ep va
        push!(line, ChU.va(bundle, exp_ξ, exp_β, :ep, rxni))
        # lb
        push!(line, ChU.lb(model, rxni))
        # ub
        push!(line, ChU.ub(model, rxni))
        # rxn eq
        push!(line, ChU.rxn_str(model, rxni))
        # rxn subSys
        push!(line, get(model.subSystems, rxni, ""))

        f(x) = replace(string(x), sep => "")
        push!(csv, join(f.(line), sep))
    end
    
    return csv
end
let exp = "E",
    csv_file = joinpath(M.MODEL_PROCESSED_DATA_DIR, "results_exp_$exp.csv")

    dat = to_csv(exp)
    open(csv_file, "w") do io
        for line in dat
            println(io, line)
        end
        println(relpath(csv_file), " created")
    end
end

## ------------------------------------------------------------------
# EP AV vs VA
function av_sv_va()
    p = plot(title = fig_title, 
        xlabel = "average", ylabel = "standar dev"
    )

    serie_fun = log10
    eps = 1e-4

    sbundles = sort(collect(bundles); by = first)
    for (exp, bundle) in sbundles
        
        exp_ξ =  Rd.val("ξ", exp)
        exp_β = closest_βs[exp]
        model = bundle[exp_ξ, :net]
        avs, vas = [], []
        for rxni in eachindex(model.rxns)
            push!(vas, ChU.va(bundle, exp_ξ, exp_β, :ep, rxni))
            push!(avs, ChU.av(bundle, exp_ξ, exp_β, :ep, rxni))
        end
        scatter!(p, serie_fun.(abs.(avs) .+ eps), 
                    serie_fun.(sqrt.(vas) .+ eps); 
                    label = "", c = colors[exp])
    end
    p
end
av_sv_va()

## ------------------------------------------------------------------
# # stoi_err_by_compartments
# # errs = log.(norm2_stoi_err(model, epout) .+ 1e-3);
# # histogram(errs; normalize = true, xlabel = "log err", ylabel = "prob dens")

# function stoi_err_by_compartments(get_comp::Function, model, out)
    
#     errs = norm2_stoi_err(model, out)

#     ## Split error by compartment
#     cerrs = Dict()
#     for (meti, met) in model.mets |> enumerate
#         c = get_comp(model, met) 
#         push!(get!(cerrs, c, []), errs[meti])
#     end

#     ps = []
#     xlims_ = [-6,4]
#     bins_ = 100
#     for (c, errs) in cerrs
#         log_err = log10.(errs .+ 1e-5) 
#         p = histogram(log_err; label = string(c), 
#             xlabel = "log err", ylabel = "prob", 
#             xlims = xlims_, bins = bins_,
#             normalize = true)
#         push!(ps, p)
#     end
#     ##
#     plot(ps...; size = [1000, 400])

# end

# cerrs = stoi_err_by_compartments(model, epout) do model, met
#     rxns = met_rxns(model, met)
#     return any(is_exchange(model, rxn) for rxn in rxns) ? "ext" : "int"
# end