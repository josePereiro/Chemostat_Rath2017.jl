#=
    This script just run all the other scripts in the correct order.
=#
#TODO: find julia native alternative for it
cd(dirname(@__FILE__)) # move to this folder
println("Now at: ", pwd())

scripts = [
    "prepare_tINIT_base_models.jl",
    "prepare_fva_pp_tINIT_models.jl",
    "fva_pp_tINIT_models_maxent_ep.jl"
]

println("To run: ", scripts)
for script in scripts
    println("\n\n----------------- Running $script -----------------\n\n")
    flush(stdout); flush(stderr)
    run(`julia --project $script`) # Change this to point to your julia executable, do not forget the flag
    flush(stdout); flush(stderr)
end