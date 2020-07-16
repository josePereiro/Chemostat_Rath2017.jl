#=
    This script just run all the other scripts in the correct order.
=#
#TODO: find julia native alternative for it
cd(dirname(@__FILE__)) # move to this folder
println("Now at: ", pwd())

# Scripts in order

base_scripts = [
    "download_raws.jl",
    "Hams_medium.jl",
    "mets_map.jl",
    "niklas_biomass.jl",
    "prepare_base_model.jl",
]

all_scripts = [
    base_scripts;
    "prepare_fva_pp_model.jl";
    "fva_pp_base_model_maxent_ep.jl";
    "fva_pp_base_model_maxent_ep_plots.jl"
]

usage = "If zero argument are passed, it will run 'all' the scripts. "*
    "Possible option 'all', 'base'."
if isempty(ARGS)
    to_run = all_scripts
elseif length(ARGS) == 1
    if ARGS[1] == "base"
        to_run = base_scripts
    elseif ARGS[1] == "all"
        to_run = all_scripts
    else
        error("$(ARGS[1]) not recognized!!\nUsage: $usage")
    end
else
    error("To many args!!\nUsage: $usage")
end

println("To run: ", to_run)
for script in to_run
    println("\n\n----------------- Running $script -----------------\n\n")
    # run(`$julia_cmd $script`)
    run(`julia --project $script`) # Change this to point to your julia executable, do not forget the flag
end