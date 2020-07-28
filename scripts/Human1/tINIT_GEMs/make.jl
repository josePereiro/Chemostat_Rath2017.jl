#TODO: find julia native alternative for it
#=
    This script just run all the other scripts in the correct order.
=#
cd(dirname(@__FILE__)) # move to this folder
println("\nNow at: ", pwd())

#######################################################
## SET UP
#######################################################
using ArgParse
set = ArgParseSettings()
@add_arg_table! set begin
    "--run", "-r", "-x"
        help = "either \"all\" (run all the scripts) or \"base\" (run only the base model scripts)"
        default = "base"
        required = false
        range_tester = (x -> x == "all" || x == "base")
    "--clear", "-c"
        help = "either \"all\" (clear all the scripts targets) or \"base\" (clear only the base model scripts targets)"
        required = false
        range_tester = (x -> x == "all" || x == "base")
end
parsed_args = parse_args(set)

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

#######################################################
## Clear
#######################################################
if !isnothing(parsed_args["clear"])
    to_clear = parsed_args["clear"] == "all" ? all_scripts : base_scripts
    println("\nTo clear: ", [basename(script.name) for script in to_clear])
    for (script, targets) in to_clear 
        for target in targets
            if isfile(target) || isdir(target)
                rm(target, force = true, recursive = true)
                println(basename(target), " deleted!!!")
            end
        end
    end
    println()
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

to_run = parsed_args["run"] == "all" ? all_scripts : base_scripts

julia = Base.julia_cmd()
println("To run: ", [basename(script.name) for script in to_run])
for (script, targets) in to_run
    println("\n\n----------------- Script $script -----------------\n")
    check_targets(targets) && continue
    flush(stdout); flush(stderr)
    run(`$julia --project $script`) # Change this to point to your julia executable, do not forget the flag
    flush(stdout); flush(stderr)
end