#=
    This script just run all the other scripts in the correct order.
=#
#TODO: find julia native alternative for it
cd(dirname(@__FILE__)) # move to this folder
println("Now at: ", pwd())

scripts = [
    "convert_rath_data.jl",
    "a1at_demand.jl"
]

println("To run: ", scripts)
for script in scripts
    println("\n\n----------------- Running $script -----------------\n\n")
    flush(stdout); flush(stderr)
    run(`julia --project $script`) # Change this to point to your julia executable, do not forget the flag
    flush(stdout); flush(stderr)
end