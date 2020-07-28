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
        help = "possible values: \"all\" (clear all the scripts targets), " *
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
clear_args = split(parsed_args["clear"], ",")
run_arg = parsed_args["run"]

cd(dirname(@__FILE__)) # move to this folder
println("\nNow at: ", pwd())

import Chemostat_Rath2017
tIG = Chemostat_Rath2017.tINIT_GEMs;

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# name-targets
targets_dict = Dict()
targets_dict["prepare_tINIT_base_models.jl"] = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, 
                                                            ["base_model_GTEx-brain-1.jls", "base_model_TCGA-GBM NT-1.jls", 
                                                            "base_model_TCGA-GBM TP-1.jls", "base_model_TCGA-GBM TR-1.jls", 
                                                            "base_model_TCGA-LGG TP-1.jls", "base_model_TCGA-LGG TR-1.jls"])
targets_dict["prepare_fva_pp_tINIT_models.jl"] = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, 
                                                            ["fva_pp_model_GTEx-brain-1.jls", "fva_pp_model_TCGA-GBM NT-1.jls", 
                                                            "fva_pp_model_TCGA-GBM TP-1.jls", "fva_pp_model_TCGA-GBM TR-1.jls", 
                                                            "fva_pp_model_TCGA-LGG TP-1.jls", "fva_pp_model_TCGA-LGG TR-1.jls"])
targets_dict["fva_pp_tINIT_models_maxent_ep.jl"] = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, 
                                                            ["fva_pp_tINIT_models_maxent_ep___GTEx-brain-1___boundles.jls",
                                                            "fva_pp_tINIT_models_maxent_ep___TCGA-GBM NT-1___boundles.jls",
                                                            "fva_pp_tINIT_models_maxent_ep___TCGA-GBM NT-1___boundles.jls",
                                                            "fva_pp_tINIT_models_maxent_ep___TCGA-GBM TR-1___boundles.jls",
                                                            "fva_pp_tINIT_models_maxent_ep___TCGA-LGG TP-1___boundles.jls",
                                                            "fva_pp_tINIT_models_maxent_ep___TCGA-LGG TR-1___boundles.jls"])
targets_dict["cache"] = joinpath.(tIG.MODEL_CACHE_DATA_DIR, readdir(tIG.MODEL_CACHE_DATA_DIR))

# Scripts-targets in order
base_scripts = ["prepare_tINIT_base_models.jl"] 
all_scripts = [base_scripts; "prepare_fva_pp_tINIT_models.jl"; "fva_pp_tINIT_models_maxent_ep.jl"]

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_args)
    to_clear = []
    for clear_arg in clear_args
        s =  clear_arg == "all" ? [all_scripts; "cache"] :
            clear_arg == "base" ? base_scripts : 
            clear_arg == "maxent_ep" ? ["fva_pp_tINIT_models_maxent_ep.jl"] :
            clear_arg == "fva_pp" ? ["prepare_fva_pp_tINIT_models.jl"] : 
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