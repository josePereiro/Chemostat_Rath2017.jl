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
        help = "possible values: \"all\" (run all the scripts),  " *
                                 "\"base\" (run only the base model scripts), " *
                                 "\"none\" (run nothing)"
        default = "base"
        required = false
        range_tester = (x -> x in ["all", "base", "none"])
    "--clear", "-c"
        help = "possible values: \"all\" (clear all the scripts targets), " *
                                "\"base\" (clear only the base model scripts targets), " *
                                "\"maxent_ep\" (clear only the maxent_ep bondles), " *
                                "\"fva_pp\" (clear only the fva preprocess models)"
        required = false
        range_tester = (x -> x in ["all", "base", "maxent_ep", "fva_pp"])
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
clear_arg = parsed_args["clear"]
run_arg = parsed_args["run"]

cd(dirname(@__FILE__)) # move to this folder
println("\nNow at: ", pwd())

import Chemostat_Rath2017
tIG = Chemostat_Rath2017.tINIT_GEMs;

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# Scripts-targets in order
base_scripts = [
    (
        name = "prepare_tINIT_base_models.jl",
        targets = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, 
            ["base_model_GTEx-brain-1.jls", "base_model_TCGA-GBM NT-1.jls", 
            "base_model_TCGA-GBM TP-1.jls", "base_model_TCGA-GBM TR-1.jls", 
            "base_model_TCGA-LGG TP-1.jls", "base_model_TCGA-LGG TR-1.jls"])
    )
] 

all_scripts = [ 
    base_scripts;
    (
        name = "prepare_fva_pp_tINIT_models.jl",
        targets = joinpath.(tIG.MODEL_PROCESSED_DATA_DIR, 
            ["fva_pp_model_GTEx-brain-1.jls", "fva_pp_model_TCGA-GBM NT-1.jls", 
            "fva_pp_model_TCGA-GBM TP-1.jls", "fva_pp_model_TCGA-GBM TR-1.jls", 
            "fva_pp_model_TCGA-LGG TP-1.jls", "fva_pp_model_TCGA-LGG TR-1.jls"])
    );
    # TODO: set this targets
    (
        name = "fva_pp_tINIT_models_maxent_ep.jl",
        targets = []
    )
]
get_names(scripts) = [basename(script.name) for script in scripts]
get_script(name) = all_scripts[get_names(all_scripts) .== name]

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_arg)

    to_clear =  clear_arg == "all" ? all_scripts :
                clear_arg == "base" ? base_scripts : 
                clear_arg == "maxent_ep" ? get_script("fva_pp_tINIT_models_maxent_ep.jl") :
                clear_arg == "fva_pp" ? get_script("prepare_fva_pp_tINIT_models.jl") : []
    
    println("\nTo clear: ", get_names(to_clear))
    for (script, targets) in to_clear 
        for target in targets
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
println("\nTo run: ", get_names(to_run))
for (script, targets) in to_run
    println("\n\n----------------- Script $script -----------------\n")
    check_targets(targets) && continue
    flush(stdout); flush(stderr)
    !dry_run_flag && run(`$julia --project $script`)
    flush(stdout); flush(stderr)
end