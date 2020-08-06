#TODO: find julia native alternative for it
#=
    This script just run all the other scripts in the correct order.
=#

#######################################################
## SET UP
#######################################################
using ArgParse

function check_clear_args(args)
    for arg in split(args, ",")
        !(arg in ["all", "base", "maxent_ep", "fva_pp", "cache"]) && return false
    end
    return true
end

set = ArgParseSettings()
@add_arg_table! set begin
    "--dry-run"
        help = "run without consecuences, just printing"
        action = :store_true
    "--run", "-r", "-x"
        help = "possible values: \"all\" (run all the scripts),  " *
                                 "\"base\" (run only the base model scripts), " *
                                 "\"none\" (run nothing)"
        default = "base"
        required = false
        range_tester = (x -> x in ["all", "base", "none"])
    "--clear", "-c"
        help = "possible values: \"raw\" (clear raw data folder), " *
                                "\"all\" (clear all the scripts targets), " *
                                "\"base\" (clear only the base model scripts targets), " *
                                "\"maxent_ep\" (clear only the maxent_ep boundles), " *
                                "\"fva_pp\" (clear only the fva preprocess models)" *
                                "\"cache\" (clear the cache forder)" *
                                "You can pass several using comma Ex: --clear=cache,maxent"
        required = false
        range_tester = check_clear_args
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
clear_args = parsed_args["clear"]
clear_args = isnothing(clear_args) ? nothing : split(clear_args, ",")
run_arg = parsed_args["run"]

using DrWatson
@quickactivate "Chemostat_Rath2017"

import Chemostat_Rath2017: HumanGEM, RAW_DATA_DIR
const HG = HumanGEM

cd(dirname(@__FILE__)) # move to this folder
println("\nNow at: ", pwd())

# name-targets
targets_dict = Dict()
targets_dict["download_raws.jl"] = joinpath.(RAW_DATA_DIR, 
                                    ["Human1_Publication_Data_Scripts", "Human1_Publication_Data_Scripts.zip"])
targets_dict["Hams_medium.jl"] = [HG.HAM_MEDIUM_FILE]
targets_dict["mets_map.jl"] = [HG.METS_MAP_FILE]
targets_dict["niklas_biomass.jl"] = [HG.NIKLAS_BIOMASS_FILE]
targets_dict["prepare_base_model.jl"] = [HG.BASE_MODEL_FILE, HG.BASE_INTAKE_INFO_FILE, 
                                    HG.EXCH_MET_MAP_FILE, HG.BASE_READABLE_MET_IDS_FILE]
targets_dict["prepare_fva_pp_model.jl"] = [HG.FVA_PP_BASE_MODEL_FILE]
targets_dict["fva_pp_base_model_maxent_ep.jl"] = []
targets_dict["fva_pp_base_model_maxent_ep_plots.jl"] = []
targets_dict["cache"] = joinpath.(HG.MODEL_CACHE_DATA_DIR, readdir(HG.MODEL_CACHE_DATA_DIR))


# scripts in order
base_scripts = ["mets_map.jl", "Hams_medium.jl", "niklas_biomass.jl", "prepare_base_model.jl"]
all_scripts = [base_scripts; "fva_pp_base_model_maxent_ep.jl"; "fva_pp_base_model_maxent_ep_plots.jl"]

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_args)
    to_clear = []
    for clear_arg in clear_args
        s =  clear_arg == "raw" ? ["download_raws.jl"] :
            clear_arg == "all" ? [all_scripts; "cache"] :
            clear_arg == "base" ? base_scripts : 
            clear_arg == "maxent_ep" ? ["fva_pp_base_model_maxent_ep.jl", "fva_pp_base_model_maxent_ep_plots.jl"] :
            clear_arg == "fva_pp" ? ["prepare_fva_pp_model.jl"] : 
            clear_arg == "cache" ? ["cache"] : []
        push!(to_clear, s...)
    end

    println("\nTo clear: ", to_clear)
    for k in to_clear 
        for target in targets_dict[k]
            if isfile(target) || isdir(target)
                !dry_run_flag && rm(target, force = true, recursive = true)
                println(basename(target), " deleted!!!")
            end
        end
    end
    flush(stdout); flush(stderr)
end

#######################################################
## RUN SCRIPTS
#######################################################

function check_targets(targets) 
    if isempty(targets)
        println("Running: No targets set\n")
        return false
    end
    for target in targets
        if !isfile(target) && !isdir(target)
            println("Running: ", basename(target), " missing\n")
            return false
        end
    end
    println("Not running: All targets present\n")
    return true
end

to_run = run_arg == "all" ? all_scripts : 
        run_arg == "base" ? base_scripts : 
        #= none =# []

julia = Base.julia_cmd()
println("\nTo run: ", to_run)
for script in to_run
    targets = targets_dict[script]
    println("\n\n----------------- Script $script -----------------\n")
    check_targets(targets) && continue
    flush(stdout); flush(stderr)
    !dry_run_flag && run(`$julia --project $script`)
    flush(stdout); flush(stderr)
end