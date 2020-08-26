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
src_file = ecG.FVA_PP_BASE_MODELS
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

