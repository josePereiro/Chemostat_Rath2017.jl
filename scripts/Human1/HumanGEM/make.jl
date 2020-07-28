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
                                "\"maxent_ep\" (clear only the maxent_ep boundles), " *
                                "\"fva_pp\" (clear only the fva preprocess models)" *
                                "\"cache\" (clear the cache forder)"
        required = false
        range_tester = (x -> x in ["all", "base", "maxent_ep", "fva_pp", "cache"])
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
clear_arg = parsed_args["clear"]
run_arg = parsed_args["run"]

import Chemostat_Rath2017
HG = Chemostat_Rath2017.HumanGEM

cd(dirname(@__FILE__)) # move to this folder
println("\nNow at: ", pwd())

# This just check that the script is run in the
# package enviroment
Chemostat_Rath2017.check_env();

# Scripts-targets in order
base_scripts = [
    (
        name = "download_raws.jl", 
        targets = joinpath.(Chemostat_Rath2017.RAW_DATA_DIR, 
            ["Human1_Publication_Data_Scripts", "Human1_Publication_Data_Scripts.zip"])
    ),
    ( 
        name = "Hams_medium.jl",
        targets = joinpath.(HG.MODEL_PROCESSED_DATA_DIR, 
            ["ham_medium.csv"])
    ),
    (
        name = "mets_map.jl",
        targets = joinpath.(HG.MODEL_PROCESSED_DATA_DIR, 
            ["mets_map.csv"])
    ),
    (
        name = "niklas_biomass.jl",
        targets = joinpath.(HG.MODEL_PROCESSED_DATA_DIR, 
            ["niklas_biomass.csv"])
    ),
    (
        name = "prepare_base_model.jl",
        targets = joinpath.(HG.MODEL_PROCESSED_DATA_DIR, 
            ["base_model.jls", "base_intake_info.csv", "exch_met_map.csv", "readable_met_ids.csv"])
    )
]

all_scripts = [ 
    base_scripts;
    (
        name = "prepare_fva_pp_model.jl",
        targets = joinpath.(HG.MODEL_PROCESSED_DATA_DIR, 
            ["fva_preprocessed_base_model.jls"])
    );
    # TODO: set this targets
    (
        name = "fva_pp_base_model_maxent_ep.jl",
        targets = []
    );
    (
        name = "fva_pp_base_model_maxent_ep_plots.jl",
        targets = []
    )
]

cache = (
    name = "cache",
    targets = joinpath.(HG.MODEL_CACHE_DATA_DIR, 
        readdir(HG.MODEL_CACHE_DATA_DIR))
)

get_names(scripts) = [basename(script.name) for script in scripts]
get_script(name) = all_scripts[get_names(all_scripts) .== name]

#######################################################
## CLEAR
#######################################################
if !isnothing(clear_arg)

    to_clear =  clear_arg == "all" ? all_scripts :
                clear_arg == "base" ? base_scripts : 
                clear_arg == "maxent_ep" ? get_script.(
                    ["fva_pp_base_model_maxent_ep.jl", "fva_pp_base_model_maxent_ep_plots"]) :
                clear_arg == "fva_pp" ? get_script("prepare_fva_pp_model.jl") : 
                clear_arg == "cache" ? get_script("cache") : []

    println("\nTo clear: ", [basename(script.name) for script in to_clear])
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