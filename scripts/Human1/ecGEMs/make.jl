#TODO: find julia native alternative for it
#=
    This script just run all the other scripts in the correct order.
=#
#######################################################
## SET UP
#######################################################
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--dry-run"
        help = "run without consecuences, just printing"
        action = :store_true
    "--run", "-r", "-x"
        help = "possible values: \"all\" (runs all scripts), " * 
                                "\"base\" (run all base scripts), " *
                                "\"none\" (runs nothing), " * 
                                "or pass the name of any scripts to just run it, " *
                                "You can pass several using comma Ex: --run=script2,script2. " *
                                "By default runs \"all\" the scripts."
        default = "all"
        required = false
    "--clear", "-c"
        help = "possible values: \"all\" (clear all the scripts targets), " *
                                "\"base\" (clear only the base model scripts targets), " *
                                "\"maxent_ep\" (clear only the maxent_ep boundles), " *
                                "\"fva_pp\" (clear only the fva preprocess models), " *
                                "\"cache\" (clear the cache forder). " * 
                                "The name of any script can be also passed. " *
                                "You can pass several using comma Ex: --clear=cache,maxent"
        required = false
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
clear_args = parsed_args["clear"]
clear_args = isnothing(clear_args) ? nothing : split(clear_args, ",")
run_args = parsed_args["run"]
run_args = isnothing(run_args) ? nothing : split(run_args, ",")

cd(@__DIR__) # move to this folder
println("\nNow at: ", pwd())

using DrWatson
quickactivate(pwd(), "Chemostat_Rath2017")

import Chemostat_Rath2017: ecGEMs
const ecG = ecGEMs

# name-targets
targets_dict = Dict()
targets_dict["prepare_ecTemplate.jl"] = [ecG.MODEL_EC_TEMPLATE_FILE]
targets_dict["generating_raw_ecModels.jl"] = [ecG.EC_BRAIN_RAW_MODELS_FILE]
targets_dict["preparing_ecModels.jl"] = []


targets_dict["cache"] = joinpath.(ecG.MODEL_CACHE_DATA_DIR, readdir(ecG.MODEL_CACHE_DATA_DIR))

# Scripts-targets in order
base_scripts = ["prepare_ecTemplate.jl", "generating_raw_ecModels.jl", "preparing_ecModels.jl"]
all_scripts = [base_scripts; []]

function dispatch_args(args)
    res = []
    for arg in args
        s =  arg == "all" ? [all_scripts; "cache"] :
            arg == "base" ? base_scripts :
            arg == "cache" ? ["cache"] : 
            arg == "none" ? [] :
            haskey(targets_dict, arg) ? [arg] : []
        push!(res, s...)
    end
    return res
end

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_args)
    to_clear = dispatch_args(clear_args)
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


to_run = dispatch_args(run_args)
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