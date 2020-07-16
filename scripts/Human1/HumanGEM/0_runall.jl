#=
    This script just run all the other scripts in the correct order.
=#
#TODO: find julia native alternative for it
cd(dirname(@__FILE__)) # move to this folder
println("Now at: ", pwd())


# Scripts in order
scripts = [
    "1_Hams_medium.jl",
    "1_mets_map.jl",
    "1_niklas_biomass.jl",
    "2_prepare_base_model.jl",
    # "3_prepare_fva_pp_model.jl",
    # "4_fva_pp_base_model_maxent_ep.jl",
    # "5_fva_pp_base_model_maxent_ep_plots.jl"
]

for script in scripts
    
    println("\n\n----------------- Running $script -----------------\n\n")
    # run(`$julia_cmd $script`)
    run(`julia --project $script`) # Change this to point to your julia executable, do not forget the flag
end