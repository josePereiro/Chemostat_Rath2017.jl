#= Model from Robinson, Jonathan L., Pınar Kocabaş, 
Hao Wang, Pierre-Etienne Cholley, Daniel Cook, Avlant Nilsson, Mihail Anton, et al. 
“An Atlas of Human Metabolism.” Science Signaling 13, no. 624 (March 24, 2020). 
https://doi.org/10.1126/scisignal.aaz1482.
=#
# Downloaded from https://github.com/SysBioChalmers/Human-GEM v1.4.0

module Human1
    import ..Chemostat_Rath2017: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR, FIGURES_DATA_DIR, RathData

    include("meta.jl")
    include("dir_and_files.jl")

end