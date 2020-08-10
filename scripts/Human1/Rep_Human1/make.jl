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
                                "\"none\" (runs nothing), " * 
                                "or pass the name of any scripts to just run it, " *
                                "You can pass several using comma Ex: --run=script2,script2. " *
                                "By default runs \"all\" the scripts."
        default = "all"
        required = false
    "--clear", "-c"
        help = "possible values: \"all\" (clear all the scripts targets), " *
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

import Chemostat_Rath2017: Rep_Human1
const RepH1 = Rep_Human1

# name-targets
targets_dict = Dict()
targets_dict["prepare_compareFVA_humanGEM_input.jl"] = [RepH1.COMP_FVA_HG_INPUT_FILE]
targets_dict["compareFVA_humanGEM.jl"] = [RepH1.COMP_FVA_HG_OUTPUT_FILE]
targets_dict["extract_data_compareFVA_humanGEM.jl"] = [RepH1.COMP_FVA_HG_EXTRACTED_DATA_FILE]

targets_dict["cache"] = joinpath.(RepH1.MODEL_CACHE_DATA_DIR, readdir(RepH1.MODEL_CACHE_DATA_DIR))

# Scripts-targets in order
all_scripts = ["prepare_compareFVA_humanGEM_input.jl", "compareFVA_humanGEM.jl",
                "extract_data_compareFVA_humanGEM.jl"]

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_args)
    to_clear = []
    for clear_arg in clear_args
        s =  clear_arg == "all" ? [all_scripts; "cache"] :
            clear_arg == "cache" ? ["cache"] : 
            haskey(targets_dict, clear_arg) ? [clear_arg] : []
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


to_run = []
for run_arg in run_args
    s =  run_arg == "all" ? [all_scripts; "cache"] : 
        run_arg == "none" ? [] : 
        haskey(targets_dict, run_arg) ? [run_arg] : []
    push!(to_run, s...)
end

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