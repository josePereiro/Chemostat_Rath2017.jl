using DrWatson

import MAT
using SparseArrays
using Test
using Distributions
using ProgressMeter

import Chemostat
import Chemostat.Utils: MetNet, rxnindex, metindex, uncompress_model
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Human1, print_action, 
                            load_cached, save_cache, set_cache_dir,
                            delete_temp_caches, temp_cache_file, RathData
import Chemostat_Rath2017.Human1: OBJ_IDER, ATPM_IDER, PROT_POOL_EXCHANGE
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;
const Rd = RathData
                            

## Loading base models
println("\nLoading ec base models")
src_file = ecG.EC_BRAIN_BASE_MODELS_FILE
@time ec_models = wload(src_file)[DATA_KEY]
println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")

## Simulation parameters
ξs = 10.0 .^ range(1, 3, length = 20)
ξs = [ξs; Rd.val(:ξ, Rd.ststs)] |> sort!
ξs_count = length(ξs);

## FBA
fbaouts_dict = Dict()
for (model_id, ec_model) in ec_models
    dict = get!(fbaouts_dict, model_id, Dict());
    model = uncompress_model(ec_model);

    # Ch.Utils.bounds!(model, PROT_POOL_EXCHANGE, 0.0, 0.072)
    for stst in Rd.ststs |> sort
        if stst == "B" || stst == "C"
            dict[stst] = dict["A"]
            continue
        end

        intake_info = HG.stst_base_intake_info(stst) |> deepcopy

        dict[stst] = []
        prog = Progress(ξs_count; desc = rpad("FBA ", 20))
        for ξ in ξs
            Ch.SteadyState.apply_bound!(model, ξ, HG.base_intake_info;
                emptyfirst = true, ignore_miss = true)
            fbaout = Ch.LP.fba(model, OBJ_IDER, PROT_POOL_EXCHANGE)
            push!(dict[stst], fbaout)
            growth = Ch.Utils.av(model, fbaout, OBJ_IDER)
            cost = Ch.Utils.av(model, fbaout, PROT_POOL_EXCHANGE)

            next!(prog; showvalues = [
                (:model, model_id),
                (:stst, stst),
                (:ξ, ξ),
                (:growth, growth),
                (:cost, cost)
            ])
        end
        finish!(prog)
    end
end

## Ploting

using Plots
pyplot();

ider_map = Dict(rath_ider => HG.exch_met_map[HG.mets_map[rath_ider]] for rath_ider in Rd.msd_mets)
ider_map["μ"] = OBJ_IDER
for (id1, id2) in ider_map
    ider_map[id2] = id1
end


## FBA vs ξ
p = plot(title = model_id, xlabel = "ξ", ylabel = "flx", xscale = :log10)
for (stst, fbaouts) in fbaouts_dict
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
scatter!(p, Rd.val(:ξ, Rd.ststs), Rd.val(:μ, Rd.ststs), label = "")
p

## Correlations
lims = [-0.1, 0.2]
p = plot(
    ylim = lims,
    xlim = lims,
    title = "FBA correlation",
    xlabel = "exp",
    ylabel = "model"
)
for (model_id, ec_model) in ec_models
    fbaouts = fbaouts_dict[model_id]
    for (stst, fbaouts) in fbaouts
        exp_ξ =  Rd.val(:ξ, stst)
        # exp_ξ
        exp_ξi = findfirst(isequal(exp_ξ), ξs)
        fbaout = fbaouts[exp_ξi]
        for rider in [Rd.msd_mets; "μ"]
            mider = ider_map[rider]
            exp_av = -Rd.qval(rider, stst)
            exp_err = Rd.qerr(rider, stst)
            mod_av = Ch.Utils.av(ec_model, fbaout, mider)
            scatter!(p, [exp_av], [mod_av]; xerr = [exp_err], label = "")
        end
    end
end
p

## Test
