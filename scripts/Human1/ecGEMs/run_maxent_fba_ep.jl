## ------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w"
        help = "number of workers to use"
        default = 1
    "--init_clear"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--finish_clear"
        help = "clear cache at the end"   
        action = :store_true
end
parsed_args = parse_args(set)
wcount = parse(Int, parsed_args["w"])
init_clear_flag = parsed_args["init_clear"]
finish_clear_flag = parsed_args["finish_clear"]

## ------------------------------------------------------------------
using Distributed

NO_WORKERS = min(length(Sys.cpu_info()) - 1, wcount)
length(workers()) < NO_WORKERS && 
    addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

## Loading everywhere
@everywhere begin

    using DrWatson 
    quickactivate(@__DIR__, "Chemostat_Rath2017")


    using Distributed
    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    import Chemostat.Utils: MetNet, EPModel,
                            rxnindex, metindex, compress_dict, 
                            uncompress_dict, clampfileds!, well_scaled_model,
                            ChstatBoundle, norm_abs_stoi_err, av, va, nzabs_range,
                            struct_to_dict

    import Chemostat.SimulationUtils: epoch_converge_ep!, cached_simulation, set_cache_dir, 
                            tagprintln_inmw, println_inmw, tagprintln_ifmw, println_ifmw,
                            save_cache, load_cache, delete_temp_caches
    import Chemostat.SteadyState: apply_bound!

    # Dev
    import Chemostat.LP: fba
    import Chemostat.Test: toy_model, simple_toy_MetNet
    import Chemostat.Test: empty_epout
    import Chemostat.MaxEntEP: maxent_ep, converge_ep!

    import Chemostat_Rath2017
    import Chemostat_Rath2017: DATA_KEY, Human1, RathData, load_data, save_data
    import Chemostat_Rath2017.Human1: OBJ_IDER, ATPM_IDER, PROT_POOL_EXCHANGE, 
                                    MAX_BOUND, ZEROTH
    const RepH1 = Human1.Rep_Human1;
    const ecG = Human1.ecGEMs
    const tIG = Human1.tINIT_GEMs;
    const HG = Human1.HumanGEM;
    const Rd = RathData
    set_cache_dir(ecG.MODEL_CACHE_DATA_DIR)
    
end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    tagprintln_inmw("CLEARING CACHE ")
    delete_temp_caches()
    println_inmw("\n")
end

## ------------------------------------------------------------------
# GLOBAL PARAMS
@everywhere begin

    const sim_params = Dict()
    sim_params[:epochlen] = 10 # This determine how often EP results will be cached

    const scaling_params = Dict()
    scaling_params[:scale_base] = 1000.0 # A smaller base could kill the process because of memory usage
    
    const epmodel_kwargs = Dict()
    epmodel_kwargs[:alpha] = 1e7 

    const epconv_kwargs = Dict()
    epconv_kwargs[:maxiter] = Int(1e4) # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:epsconv] = 1e-5 # The error threshold of convergence
    epconv_kwargs[:maxvar] = 1e10
    epconv_kwargs[:minvar] = 1e-10

    params_hash = hash((sim_params, scaling_params, epmodel_kwargs, epconv_kwargs))

end

## ------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "MAXENT_FBA_EP_v2"

## ------------------------------------------------------------------
# LOAD MODELS
tagprintln_inmw("LOADING EC MODELS")
src_file = ecG.FVA_PP_BASE_MODELS
ec_models = load_data(src_file)
model_ids = ec_models |> keys |> collect
for (model_id, model_dict) in ec_models
    local model = model_dict |> compress_dict |> MetNet
    ec_models[model_id] = model
    clampfileds!(model, [:lb, :ub, :b]; abs_max = MAX_BOUND, zeroth =  ZEROTH)
    println_ifmw("model: ", model_id, " size: ", size(model), " S nzabs_range: ", nzabs_range(model.S), "\n")
end    

## ------------------------------------------------------------------
# SCALE MODELS
# tagprintln_inmw("SCALING MODELS")
# for (model_id, model) in ec_models
#     ec_models[model_id] = well_scaled_model(model, scaling_params[:scale_base]; verbose = false)
#     println_ifmw("model: ", model_id, " size: ", size(model), " S nzabs_range: ", nzabs_range(model.S), "\n")
# end  

## ------------------------------------------------------------------
# CACHE MODELS
@everywhere models_cache_id = (:MODELS, sim_global_id)
save_cache(models_cache_id, ec_models; headline = "MODELS CACHE SAVED")
# free 
ec_models = nothing
GC.gc()

## ------------------------------------------------------------------
# GET MODEL FUNCTION
@everywhere function prepare_model(model_id, ξ, stst)
    dat = load_cache(models_cache_id; verbose = false)
    isnothing(dat) && error("Unable to load model!!")
    model = dat[model_id] 

    # intake info
    intake_info = HG.stst_base_intake_info(stst)

    # Chemostat steady state constraint, see Cossio's paper, (see README)
    apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)

    return model
end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
# delete_temp_caches()

## ------------------------------------------------------------------
# SIMULATION
# Any of the loops can be parallelized by just changing one of the 'map' functions

# This uniquely identify a worker
@everywhere indexid(wid) = (:INDEX, wid, sim_global_id)

to_map = Iterators.product(model_ids)
map(model_ids) do (model_id)

    tagprintln_inmw("PROCESSING MODEL ", 
        "\nid: ", model_id, 
        "\n")
    
    to_map = Iterators.product(Rd.ststs, [model_id])
    map(to_map) do (stst, model_id)
        
        ## SIMULATION PARAMS
        ξs = [Rd.val(:ξ, stst)]
        # ξs = rand(5) # Test
        βs = [0.0; range(1e3, 1e5, length = 25)] 
        # βs = [0.0] # Test
        
        to_map = Iterators.product(ξs, [βs], [stst], [model_id])
        pmap(to_map) do (ξ, βs, stst, model_id)
            
            ## HASH SEEDS
            # TODO: add PARAMS hash
            model_hash = (model_id, sim_global_id)
            stst_hash = (stst, model_hash)
            sim_hash = (ξ, stst_hash)

            dat = cached_simulation(;
                epochlen = sim_params[:epochlen], 
                # epochlen = 100, # Test
                verbose = true,
                sim_id = sim_hash,
                get_model = function()
                    return prepare_model(model_id, ξ, stst);
                    # return simple_toy_MetNet() # Test
                end,
                objider = OBJ_IDER, 
                # objider = "biom", # Test
                beta_info = [(OBJ_IDER, βs)],
                # beta_info = [("biom", βs)], 
                costider = PROT_POOL_EXCHANGE,
                clear_cache = false,
                use_seed = true,
                epmodel_kwargs = epmodel_kwargs,
                epconv_kwargs = epconv_kwargs
            )
            
            ## SAVING DATA
            # model = prepare_model(model_id, ξ, stst)
            model = simple_toy_MetNet()
            save_cache(sim_hash, (model_id, stst, ξ, βs, model, dat); 
                headline = "CATCHING RESULTS\n")
            
            ## SAVING TO INDEX
            # Here I save a link to the work done by worker
            index_id = indexid(myid())
            index = load_cache(index_id; verbose = false)
            index = isnothing(index) ? [] : index
            push!(index, sim_hash)
            save_cache(index_id, index; headline = "UPDATING INDEX")
            
            GC.gc()
            return nothing
        end # map(ξs) do ξ

        return nothing
    end # map(Rd.ststs) do stst

    return nothing
end # map(model_ids) do model_id

## COLLECTING RESULTS
tagprintln_inmw("COLLECTING RESULTS ")
boundles = Dict()
for wid in workers()

    index_id = indexid(wid)
    index = load_cache(index_id; verbose = false)
    isnothing(index) && continue
    
    for id in index
        
        model_id, stst, ξ, βs, model, dat = load_cache(id; verbose = false)
       
        # Bundle
        model_dict = get!(boundles, model_id, Dict())
        boundle = get!(model_dict, stst, ChstatBoundle())

        boundle[ξ, :net] = model
        boundle[ξ, :fba] = dat[:fba]

        for (βi, β) in βs |> enumerate
            boundle[ξ, β, :ep] = dat[(:ep, βi)]
        end

    end
end

## SAVING
tagprintln_inmw("SAVING RESULTS ")
save_data(ecG.MAXENT_FBA_EB_BOUNDLES_FILE, boundles)

## CLEAR CACHE (WARNING)
if finish_clear_flag
    tagprintln_inmw("CLEARING CACHE ")
    delete_temp_caches()
end