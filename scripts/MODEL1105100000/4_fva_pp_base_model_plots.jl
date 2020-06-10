# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

using Revise
using DataFrames
using Serialization
using StatsBase
using Statistics
using JLD
using Plots
pyplot()

import Chemostat
Ch = Chemostat;

using Revise
import Chemostat_Rath2017
Rd = Chemostat_Rath2017.RathData
M = Chemostat_Rath2017.MODEL1105100000;

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# ---
# ## Meta
# ---

notebook_name = "fva_pp_base_model_plot";
data_notebook_name = "fva_pp_base_model_maxent_ep_v1";

# ---
# ## MaxEnt-EP
# --- 

# Check cached objects
cache_file = joinpath(M.MODEL_CACHE_DATA_DIR, "$(data_notebook_name)___boundles.jls");
# cache_file = "/Volumes/Store/Temp Results/Chemostat_Rath2017/fva_pp_base_model_maxent_ep_v1___boundles.jls"
params, data = deserialize(cache_file)
println("loaded $(relpath(cache_file))!!!")

boundles = Dict()
for (stst, boundle) in data
    boundles[stst] = boundle
end
boundles["B"] = boundles["A"]
boundles["C"] = boundles["A"];



# +
# ### Cleaning boundles, keep only the complete errorless indexes
# function get_clean(boundle)
#     βs_ = sort(boundle.βs, rev = true)
#     ξs_ = sort(boundle.ξs)
    
#     errorless_βi = trues(length(βs_))
#     for ξ in ξs_
#         for (βi, β) in enumerate(βs_)
#             if haskey(boundle, ξ, β, :ep)
#                 data = Ch.Utils.get_data(boundle, ξ, β, :ep)
#                 errorless_βi[βi] = data isa Ch.Utils.EPout && errorless_βi[βi]
#             else
#                 errorless_βi[βi] = false
#             end
#         end
#     end
#     all(errorless_βi .== false) && error("Boundle is too dirty!!!")
    
#     errorless_βi = βs_[errorless_βi]
#     errorless_boundle = Ch.Utils.ChstatBoundle()
#     for ξ in ξs_
#         model = Ch.Utils.get_data(boundle, ξ, :net)
#         fbaout = Ch.Utils.get_data(boundle, ξ, :fba)
#         Ch.Utils.add_data!(errorless_boundle, ξ, :net, model)
#         Ch.Utils.add_data!(errorless_boundle, ξ, :fba, fbaout)
        
#         for β in errorless_βi
#             epout = Ch.Utils.get_data(boundle, ξ, β, :ep)
#             Ch.Utils.add_data!(errorless_boundle, ξ, β, :ep, epout)
#         end
#     end
#     return errorless_boundle
# end

# +
# for (stst, boundle) in boundles
#     boundles[stst] = get_clean(boundle)
# end
# # This three ststs are equivalents (repetitions of the same experiment)
# # boundles["A"] = boundles["C"]
# # boundles["B"] = boundles["C"];
# -

obj_ider = params["obj_ider"]; # all models are equals in this sense

# ### Closest βs
# We select the MaxEnt distribution that better explaint the experimental mu

### maxent βi
closest_βs = Dict()
for (stst, boundle) in boundles

    exp_μ = Rd.val("μ", stst)
    exp_ξ = Rd.val("ξ", stst)

    closest_β = boundle.βs[1]
    for β in boundle.βs
        curr_μ = Ch.Utils.av(boundle, exp_ξ, β, :ep, obj_ider)
        closest_μ = Ch.Utils.av(boundle, exp_ξ, closest_β, :ep, obj_ider)
        closest_β = abs(curr_μ - exp_μ) < abs(closest_μ - exp_μ) ? β : closest_β
    end
    
    closest_ep_μ = Ch.Utils.av(boundle, exp_ξ, closest_β, :ep, obj_ider)
    println("stst: $stst")
    println("exp_mu: $exp_μ,  closest_mu: $closest_ep_μ  rel exp/model: $(exp_μ/closest_ep_μ)")
    closest_βs[stst] = closest_β
end

# ---
# ## Plots
# ---

exchs = [M.exch_met_map[M.mets_map[met]] for met in Rd.msd_mets];

colors = Dict()
for (stst, c) in zip(Rd.ststs, Plots.distinguishable_colors(length(Rd.ststs)))
    colors[stst] = c
end
colors["A"] = colors["C"];
colors["B"] = colors["C"];

# ### Plot flxs as function of xi

# +
function plot_xis!(p, ider, ststs::Vector = collect(keys(boundles)))
    for stst in ststs
        boundle = boundles[stst]
        sense = ider == "BIOMASS" ? 1 : -1

        # fba
        fba_avs = Ch.Utils.av(boundle, boundle.ξs, :fba, ider)
        Plots.plot!(p, boundle.ξs, sense .* fba_avs, color = colors[stst], label = "", ls = :dash)

        # ep
        β = closest_βs[stst]
        βstr = floor(Int, β)
        ep_avs = Ch.Utils.av(boundle, boundle.ξs, β, :ep, ider)
        ep_stds = sqrt.(Ch.Utils.va(boundle, boundle.ξs, β, :ep, ider))
        Plots.plot!(p, boundle.ξs, sense .* ep_avs, 
            yerr = ep_stds,
            color = colors[stst], label = "")

    end
    
    for stst in Rd.ststs
        # exp
        exp_ider = ider == "BIOMASS" ? "μ" : "q" * M.mets_map[M.exch_met_map[ider]]


        exp_av = Rd.val(exp_ider, stst)
        exp_av_err = Rd.err(exp_ider, stst)

        exp_ξ = Rd.val("ξ", stst)
        exp_ξ_err = Rd.err("ξ", stst)

        Plots.scatter!(p, [exp_ξ], [exp_av], xerr = [exp_ξ_err], yerr = [exp_av_err], 
            color = colors[stst], ms = 10, label = "")
    end
    return p
end

function plot_xis_legend(ststs = Rd.ststs; kwargs...)
    p = Plots.plot(;legend = :topleft, framestyle = :none, kwargs...)
    
    Plots.plot!(p, [], [], ls = :dash, color = :black, label = "fba")
    Plots.plot!(p, [], [], color = :black, label = "ep")
    
    # exp colors
    for (stst, color) in colors
        Plots.scatter!(p, [], [], m = (15, :circle), label = "exp_$stst", color = color)
    end
    return p
end
# -

ps = []
for ider in [obj_ider; exchs]
    p = Plots.plot(title = ider,  xlabel = "xi",  ylabel = "flx", legend = false)
    push!(ps, plot_xis!(p, ider))
end
push!(ps, plot_xis_legend());

p = Plots.plot(ps..., size = [1200, 900], titlefont = 9) 

# +
# fig_file = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(data_notebook_name)___flx_vs_xi.png")
# Plots.savefig(p, fig_file)
# println(relpath(fig_file), " created!!!")
# -

# ### Cost

# +
ider = "enzyme_solvent_capacity"
p = Plots.plot(title = ider,  xlabel = "xi",  ylabel = "flx", legend = false)
for stst in Rd.ststs
    boundle = boundles[stst]

    # fba
    fba_avs = Ch.Utils.av(boundle, boundle.ξs, :fba, ider)
    Plots.plot!(p, boundle.ξs, fba_avs, color = colors[stst], label = "", ls = :dash)

    # ep
    β = closest_βs[stst]
    βstr = floor(Int, β)
    ep_avs = Ch.Utils.av(boundle, boundle.ξs, β, :ep, ider)
    ep_stds = sqrt.(Ch.Utils.va(boundle, boundle.ξs, β, :ep, ider))
    Plots.plot!(p, boundle.ξs, ep_avs, 
        yerr = ep_stds,
        color = colors[stst], label = "")

end
p
# -

# ### Total corr

# +
p = Plots.plot()
function plot_tot_corr!(fba_p, ep_p, iders = [obj_ider; exchs], 
        ststs = collect(keys(boundles)))
    
    xlim_ = (lb = Inf, ub = -Inf)
    ylim_ = (lb = Inf, ub = -Inf)
    
    for ider in iders
        for stst in ststs
            boundle = boundles[stst]

            # exp
            exp_ider = ider == "BIOMASS" ? "μ" : "q" * M.mets_map[M.exch_met_map[ider]]
            sense = ider == "BIOMASS" ? 1 : -1

            exp_av = sense * Rd.val(exp_ider, stst)
            exp_av_err = Rd.err(exp_ider, stst)

            exp_ξ = Rd.val("ξ", stst)


            # fba
            fba_av = Ch.Utils.av(boundle, exp_ξ, :fba, ider)
            Plots.scatter!(fba_p, [exp_av], [fba_av], xerr = exp_av_err,  
                m = (3, :square), label = "", color = colors[stst])

#             # ep
            closest_β = closest_βs[stst]
            ep_av = Ch.Utils.av(boundle, exp_ξ, closest_β, :ep, ider)
            ep_std = sqrt(Ch.Utils.va(boundle, exp_ξ, closest_β, :ep, ider))
            Plots.scatter!(ep_p, [exp_av], [ep_av], xerr = exp_av_err, yerr = ep_std, 
                m = (3, :star), label = "", color = colors[stst])

#             println("stst: $stst")
#             println("exp_mu: $exp_av,  closest_mu: $ep_av  rel exp/model: $(exp_av/ep_av)")

#             limits
            xlim_ = (lb = min(xlim_.lb, exp_av), ub = max(xlim_.ub, exp_av))
            ylim_ = (lb = minimum([ylim_.lb; fba_av; ep_av ]), 
                     ub = maximum([ylim_.ub; fba_av; ep_av]))
        end
    end
    # ref
#     Plots.plot!(ep_p, x -> x, lw = 2, ls = :dash, color = :black, label = "")
#     Plots.plot!(fba_p, x -> x, lw = 2, ls = :dash, color = :black, label = "")

#     limits
#     xmargin_ = 0.15*(xlim_.ub - xlim_.lb)
#     xlim_ = (xlim_.lb - xmargin_, xlim_.ub + xmargin_)
#     ymargin_ = 0.15*(ylim_.ub - ylim_.lb)
#     ylim_ = (ylim_.lb - ymargin_, ylim_.ub + ymargin_)
#     Plots.plot!(ep_p, xlim = xlim_, ylim = ylim_)
#     Plots.plot!(fba_p, xlim = xlim_, ylim = ylim_)
    
end
# -

fba_p = Plots.plot(title = "FBA", xlabel = "experimental", ylabel = "modeled", legend = false)
ep_p = Plots.plot(title = "EP", xlabel = "experimental", ylabel = "modeled", legend = false)
# Plots.plot!(fba_p, x -> x)
# Plots.plot!(ep_p, x -> x)
plot_tot_corr!(fba_p, ep_p)
tot_corr_p = Plots.plot(fba_p, ep_p, size = [600, 300], titlefont = 8, 
    yaxis = [-0.2, 0.2],
    xaxis = [-0.2, 0.2],
)

# +
# fig_file = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(data_notebook_name)___tot_corrs.png")
# Plots.savefig(tot_corr_p, fig_file)
# println(relpath(fig_file), " created!!!")
# -

# ### Correlations

# +
function plot_corr!(p, ider, 
        ststs = collect(keys(boundles)); ref = true)
    
    xlim_ = (lb = Inf, ub = -Inf)
    ylim_ = (lb = Inf, ub = -Inf)

    for stst in ststs
        boundle = boundles[stst]

        # exp
        exp_ider = ider == "BIOMASS" ? "μ" : "q" * M.mets_map[M.exch_met_map[ider]]
        sense = ider == "BIOMASS" ? 1 : -1
        
        exp_av = sense * Rd.val(exp_ider, stst)
        exp_av_err = Rd.err(exp_ider, stst)

        exp_ξ = Rd.val("ξ", stst)


        # fba
        fba_av = Ch.Utils.av(boundle, exp_ξ, :fba, ider)
        Plots.scatter!(p, [exp_av], [fba_av], xerr = exp_av_err,  
            m = (10, :square), label = "", color = colors[stst])

        # ep
        closest_β = closest_βs[stst]
        ep_av = Ch.Utils.av(boundle, exp_ξ, closest_β, :ep, ider)
        ep_std = sqrt(Ch.Utils.va(boundle, exp_ξ, closest_β, :ep, ider))
        Plots.scatter!(p, [exp_av], [ep_av], xerr = exp_av_err, yerr = ep_std, 
            m = (15, :star), label = "", color = colors[stst])
        
#         println("stst: $stst")
#         println("exp_mu: $exp_av,  closest_mu: $ep_av  rel exp/model: $(exp_av/ep_av)")

        # limits
        xlim_ = (lb = min(xlim_.lb, exp_av), ub = max(xlim_.ub, exp_av))
        ylim_ = (lb = minimum([ylim_.lb; fba_av; ep_av]), ub = maximum([ylim_.ub; fba_av; ep_av]))
    end
    # ref
    ref && Plots.plot!(p, x -> x, -1e3, 1e3, lw = 2, ls = :dash, color = :black, label = "")

    # limits
    xmargin_ = 0.15*(xlim_.ub - xlim_.lb)
    xlim_ = (xlim_.lb - xmargin_, xlim_.ub + xmargin_)
    ymargin_ = 0.15*(ylim_.ub - ylim_.lb)
    ylim_ = (ylim_.lb - ymargin_, ylim_.ub + ymargin_)
    Plots.plot!(xlim = xlim_, ylim = ylim_)
    
end

function plot_corr_legend(ststs = Rd.ststs; kwargs...)
    p = Plots.plot(;legend = :topleft, framestyle = :none, kwargs...)
    # data legend
    Plots.scatter!(p, [], [], m = (15, :square), label = "fba", color = :black)
    Plots.scatter!(p, [], [], m = (15, :star), label = "ep", color = :black)
    
    # exp colors
    for (stst, color) in colors
        Plots.scatter!(p, [], [], m = (15, :circle), label = "exp_$stst", color = color)
    end
    return p
end
# -

ps = []
for ider in [obj_ider; exchs]
    p = Plots.plot(title = ider, xlabel = "experimental", ylabel = "modeled", legend = false)
    push!(ps, plot_corr!(p, ider))
end
push!(ps, plot_corr_legend());

Plots.plot(ps..., size = [1200, 900], titlefont = 10)

# +
# fig_file = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(data_notebook_name)___individual_corrs.png")
# Plots.savefig(tot_corr_p, fig_file)
# println(relpath(fig_file), " created!!!")
# -

# ### Marginals

function plot_marginal!(p, ider, stst)
    ξ_ = Rd.val("ξ", stst)
    β_ = closest_βs[stst]
    
    boundle = boundles[stst]

    # marginals
    Ch.Plots.plot_marginal!(p, boundle, ξ_, β_, [:fba, :ep], ider)

    # ep mean
    ep_av = Ch.Utils.av(boundle, ξ_, β_, :ep, ider)
    ep_std = sqrt(Ch.Utils.va(boundle, ξ_, β_, :ep, ider))
    Plots.vline!(p,[ep_av], color = :red, lw = 5)
    # Plots.vline!(p,[ep_av - ep_std], color = :red, lw = 2)
    # Plots.vline!(p,[ep_av + ep_std], color = :red, lw = 2)

    # exp value
    exp_ider = ider == "BIOMASS" ? "μ" : "q" * M.mets_map[M.exch_met_map[ider]]
    sense = ider == "BIOMASS" ? 1.0 : -1.0
    exp_av = sense * Rd.val(exp_ider, stst)
    exp_av_err = Rd.err(exp_ider, stst)
    
    exp_color = :gold
    Plots.vline!(p,[exp_av], color = exp_color, lw = 3, ls = :dash)
    Plots.vline!(p,[exp_av - exp_av_err], color = exp_color, lw = 1, ls = :dash)
    Plots.vline!(p,[exp_av + exp_av_err], color = exp_color, lw = 1, ls = :dash)
    
    return p
end

# +
xlims_ = Dict()
xlims_[obj_ider] = [-0.01, 0.05]
xlims_["EX_glc_LPAREN_e_RPAREN"] = [-0.2, 0.2]
xlims_["EX_lac_L_LPAREN_e_RPAREN_"] = [-1, 10]
xlims_["EX_gln_L_LPAREN_e_RPAREN_"] = [-0.1, 0.1]
xlims_["EX_nh4_LPAREN_e_RPAREN"] = [-1, 4]
xlims_["EX_gal_LPAREN_e_RPAREN_"] = [-1, 4]
xlims_["EX_pyr_LPAREN_e_RPAREN_"] = [-1, 8]
xlims_["EX_glu_L_LPAREN_e_RPAREN_"] = [-0.1, 0.4]
xlims_["EX_ala_L_LPAREN_e_RPAREN_"] = [-0.2, 1]
# xlims_["EX_asp_L_LPAREN_e_RPAREN_"] = [-0.2, 1];

psr = []
for stst in boundles |> keys |> collect
    psc = []
    for ider in [obj_ider; exchs]
        occursin("EX_gal", ider) && continue
        p = Plots.plot(title = "$ider\nexp: $stst", 
            xlabel = "flx", ylabel = "prob", 
            xlims = get(xlims_, ider, false),
            legend = false)
        plot_marginal!(p, ider, stst)
        push!(psc, p)
    end
    p = Plots.plot(psc..., layout = grid(length(psc), 1))
    push!(psr, p)
end
# -

Plots.plot(psr..., size = [1300,2000], titlefont = 9, layout = grid(1, length(psr)))

# +
# fig_file = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(data_notebook_name)___individual_marginals.png")
# Plots.savefig(tot_corr_p, fig_file)
# println(relpath(fig_file), " created!!!")
# -

# ### Stoi err

# +
p = Plots.plot(xlabel = "beta", 
    ylabel = "abs_stoi_err/ abs_mean_flx", 
    legend = :topleft)

for (stst, boundle) in boundles
    
    ξ = Rd.val(:ξ, stst)

    metnet = boundle[ξ, :net]
    
    M,B = length(metnet.mets), length(boundle.βs)
    max_abs_errs_ = []
    mean_abs_errs_ = []
    std_abs_errs_ = []
    
    
    for β in boundle.βs 
        abs_errs_ = abs.(Ch.Utils.stoi_err(boundle, ξ, β, :ep))
        
        max_abs_err_ = maximum(abs_errs_)
        mean_abs_err_ = mean(abs_errs_)
        std_abs_err_ = std(abs_errs_)
        
        
        mean_abs_flxs_ = mean(abs.(Ch.Utils.av(boundle, ξ, β, :ep)))
        
        push!(max_abs_errs_, max_abs_err_/ mean_abs_flxs_)
        push!(mean_abs_errs_, mean_abs_err_/ mean_abs_flxs_)
        push!(std_abs_errs_, std_abs_err_/ mean_abs_flxs_)
    end
    
    Plots.plot!(p, boundle.βs , max_abs_errs_, color = colors[stst], 
        lw = 1, label = "")
    Plots.plot!(p, boundle.βs , mean_abs_errs_, yerr = std_abs_errs_, alpha = 0.5, color = colors[stst], 
        lw = 1, label = "", ls = :dash)
    
end
Plots.plot!(p, [], [], lw = 3, label = "max err", color = :black)
Plots.plot!(p, [], [], lw = 3, label = "mean/std err", ls = :dash, color = :black)
p

# +
# fig_file = joinpath(M.MODEL_FIGURES_DATA_DIR, "$(data_notebook_name)___norm_stoi_err.png")
# Plots.savefig(tot_corr_p, fig_file)
# println(relpath(fig_file), " created!!!")
# -

boundle = boundles["A"]
model = boundle[1, :net]
e = Ch.Utils.norm_abs_stoi_err(boundle, 1, 1, :ep)
for (met, e) in zip(model.mets, e)
    e < 0.1 && continue
    println(met, ", e: ", e)
#     println(Ch.Utils.balance_str(boundle, 1, 1, :ep, met, digits = 3))
#     println()
end

Ch.Utils.norm_abs_stoi_err(boundle, 1, 1, :ep, "enzyme_c")


