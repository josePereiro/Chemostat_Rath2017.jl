import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

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
    "--not-run", "-n"
        help = "run any script"
        action = :store_true
    "--clear", "-c"
        help = "Clear first all the scripts target files"
        action = :store_true
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
clear_arg = parsed_args["clear"]
not_run_arg = parsed_args["not-run"]

cd(@__DIR__) # move to this folder
println("\nNow at: ", pwd())

import Chemostat_Rath2017
const Rd = Chemostat_Rath2017.RathData;

# Scripts-targets in order
all_scripts = [
    (
        name = "convert_rath_data.jl",
        targets = joinpath.(Rd.RATH_PROCESSED_DATA_DIR, 
            ["rath2017___42_MAX_UB_standard_medium.tsv", "rath2017___cont_exp_A.tsv",
            "rath2017___cont_exp_B.tsv", "rath2017___cont_exp_C.tsv", 
            "rath2017___cont_exp_D.tsv", "rath2017___cont_exp_E.tsv",
            "rath2017___cont_exp_F01.tsv", "rath2017___cont_exp_F02.tsv",
            "rath2017___cont_exp_F03.tsv", "rath2017___max_invitro_fluxs.tsv"])
    ),
    (
        name = "a1at_demand.jl",
        targets = joinpath.(Rd.RATH_PROCESSED_DATA_DIR, 
            ["a1at_aa_rel_abundance.csv"])
    )
] 

get_names(scripts) = [basename(script.name) for script in scripts]
get_script(name) = all_scripts[get_names(all_scripts) .== name]

#######################################################
## CLEAR
#######################################################
if clear_arg

    to_clear = all_scripts
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

if !not_run_arg
    to_run = all_scripts
    julia = Base.julia_cmd()
    println("\nTo run: ", get_names(to_run))
    for (script, targets) in to_run
        println("\n\n----------------- Script $script -----------------\n")
        check_targets(targets) && continue
        flush(stdout); flush(stderr)
        !dry_run_flag && run(`$julia --project $script`)
        flush(stdout); flush(stderr)
    end
end