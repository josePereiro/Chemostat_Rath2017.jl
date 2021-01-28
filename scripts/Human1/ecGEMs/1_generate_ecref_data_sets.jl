## ----------------------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

@time begin
    import MAT
    using SparseArrays
    using Test
    using Distributions

    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const H1 = ChR.Human1
    const ecG = H1.ecGEMs
    const HG = H1.HumanGEM
end

## ----------------------------------------------------------------------------
# Collect files
# Pick all ec Models
data_dir = ecG.EC_RAW_MODELS_DIR
ec_model_files = Dict()
for (root, _, files) in walkdir(data_dir)
    for file in files
        if basename(file) == "ecModel_batch.mat"
            model_id = basename(root)
    
            srcfile = joinpath(root, model_id * ".mat")
            !isfile(srcfile) && error(relpath(srcfile), " not found")
            ecfile = joinpath(root, "ecModel_batch.mat")
            !isfile(ecfile) && error(relpath(ecfile), " not found")
            
            ec_model_files[model_id] = (src = srcfile, ec = ecfile)
            break # Test
        end
    end
end

## ----------------------------------------------------------------------------
# Processing

println("\n\n--------------- Preprocessing models ---------------")
ec_reference_models = Dict()
models_count = length(ec_model_files)
for (i, (id, files)) in ec_model_files |> enumerate

    println("\n\nmodel $id [$i/$models_count] ---------------\n")
    ec_reference_models[id] = Dict()

    src_exchs = nothing # this need to survive both iterations
    for model_sym in [:src, :ec]
        is_src = model_sym == :src

        println("\nProcessing ", model_sym)
        model = Ch.Utils.read_mat(files[model_sym]); 
        Ch.Utils.clamp_bounds!(model, Human1.MAX_BOUND, Human1.ZEROTH);
        
        # We delete the boundary metabolites, they are not required
        # bkwd and fwd splitted reactions are troublemakers for EP, but they 
        # are necessary to model enzymatic costs. So, we leave as least as possible 
        # only leaving in rev format the exchanges

        # We unified the exchanges (make them a unique rxn), and let them in a 
        # semi-open state (intake bloked, outtake open)   
        model = Human1.delete_boundary_mets(model);
        
        # We close all reactions marked as "Exchanges/boundary" and returns the exchanges
        # We also close complitly the reactions classified as exchanges by having anly reacts or
        # prods
        exchs = Human1.prepare_extract_exchanges!(model)
        is_src && (src_exchs = exchs)
        
        # We delate any backward defined exchange reaction, 
        # all reactions will be possible reversible
        !is_src && (model = Human1.del_REV_rxns(model, src_exchs))

        # Apply Hams medium
        medium = HG.base_intake_info |> keys |> collect
        medium = filter((rxn) -> rxn in model.rxns, medium)
        Human1.open_rxns!(model, medium)

        # Open prot pool exchange
        !is_src && Ch.Utils.bounds!(model, Human1.PROT_POOL_EXCHANGE, 0.0, Human1.MAX_BOUND);

        # Check that the model is feasible
        fbaout = Human1.try_fba(model, Human1.BIOMASS_IDER);
        @assert fbaout.obj_val > Human1.ZEROTH

        model = Ch.Utils.compress_model(model);
        ec_reference_models[id][model_sym] = model
    end
end

# ---
# ## Generate ecmaps

refs = Dict()
println("\n\n--------------- Generating ec reference data ---------------")
models_count = length(ec_reference_models)
for (i, (model_id, models_dict)) in ec_reference_models |> enumerate
    
    println("\n\n model $model_id [$i/$models_count]\n")
    
    src_model = models_dict[:src]
    ec_model = models_dict[:ec]

    ec_refdata = Human1.get_ec_refdata(src_model, ec_model);
    new_ec_model = Human1.build_ecModel(src_model, [ec_refdata]);
    Human1.print_ec_stats(new_ec_model)
    
    ec_model = Ch.Utils.compress_model(ec_model)
    new_ec_model = Ch.Utils.compress_model(new_ec_model)
    
    # testing
    # We used a tINIT GEM and its respective ecModel as ec template.
    # So, the resulting new ecModel must be equal to the template one
    @assert all(new_ec_model.S[:] |> sort .== ec_model.S[:] |> sort)
    @assert all(new_ec_model.b[:] |> sort .== ec_model.b[:] |> sort)
    @assert all(new_ec_model.lb[:] |> sort .== ec_model.lb[:] |> sort)
    @assert all(new_ec_model.ub[:] |> sort .== ec_model.ub[:] |> sort)
    @assert all(new_ec_model.mets[:] |> sort .== ec_model.mets[:] |> sort)
    @assert all(new_ec_model.rxns[:] |> sort .== ec_model.rxns[:] |> sort)
    println("All tests passed!!!"); flush(stderr)
    
    refs[model_id] = ec_refdata
end

# saving
file = ecG.EC_REFERENCE_DATA
tagsave(file, Dict(DATA_KEY => refs))
println(relpath(file), " created!!!, size: ", filesize(file), " bytes")


