## ----------------------------------------------------------------------------
import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_Rath2017")

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

    import UtilsJL
    const UJL = UtilsJL
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
let
    models_count = length(ec_model_files)
    base_intake_info = HG.load_base_intake_info()

    for (i, (id, files)) in ec_model_files |> enumerate

        println("\n\nmodel $id [$i/$models_count] ---------------\n")
        ec_reference_models[id] = Dict()

        src_exchs = nothing # this need to survive both iterations
        for model_sym in [:src, :ec]
            is_src = model_sym == :src

            println("\nProcessing ", model_sym)
            model = ChU.read_mat(files[model_sym])
            model = ChU.fix_dims(model)
            ChU.clampfields!(model, [:lb, :ub]; 
                abs_max = H1.MAX_BOUND, zeroth = H1.ZEROTH
            )

            
            # We delete the boundary metabolites, they are not required
            # bkwd and fwd splitted reactions are troublemakers for EP, but they 
            # are necessary to model enzymatic costs. So, we leave as least as possible 
            # only leaving in rev format the exchanges

            # We unified the exchanges (make them a unique rxn), and let them in a 
            # semi-open state (intake bloked, outtake open)   
            # global m = model
            model = H1.delete_boundary_mets(model);
            
            # We close all reactions marked as "Exchanges/boundary" and returns the exchanges
            # We also close completely the reactions classified as exchanges by having any reacts or
            # prods
            exchs = H1.prepare_extract_exchanges!(model)
            is_src && (src_exchs = exchs)
            
            # We delate any backward defined exchange reaction, 
            # all reactions will be possible reversible
            !is_src && (model = H1.del_REV_rxns(model, src_exchs))

            # Apply Hams medium
            medium = base_intake_info |> keys |> collect
            medium = filter((rxn) -> rxn in model.rxns, medium)
            H1.open_rxns!(model, medium)

            # Open prot pool exchange
            !is_src && ChU.bounds!(model, H1.PROT_POOL_EXCHANGE, 0.0, H1.MAX_BOUND);

            # Check that the model is feasible
            fbaout = H1.try_fba(model, H1.BIOMASS_IDER);
            # @assert fbaout.obj_val > H1.ZEROTH

            model = ChU.compressed_model(model);
            ec_reference_models[id][model_sym] = model
        end # for model_sym 
        break # Test
    end # for (i, (id, files)) 
end

# ## ----------------------------------------------------------------------------
# let
#     m0 = deppcopy(m)
#     # ChU.summary(m, fbaout)
#     exchs = ChU.exchanges(m0)
#     for exch in exchs
#         m0.lb[exch] == 0.0 && continue
#         # ChU.summary(m0, m0.rxns[exch])
#         ChU.
#     end
#     H1.try_fba(m0, H1.BIOMASS_IDER);
# end

## ----------------------------------------------------------------------------
# Generate ecmaps

refs = Dict()
println("\n\n--------------- Generating ec reference data ---------------")
models_count = length(ec_reference_models)
for (i, (model_id, models_dict)) in ec_reference_models |> enumerate
    
    println("\n\n model $model_id [$i/$models_count]\n")
    
    src_model = models_dict[:src]
    ec_model = models_dict[:ec]

    ec_refdata = H1.get_ec_refdata(src_model, ec_model);
    new_ec_model = H1.build_ecModel(src_model, [ec_refdata]);
    H1.print_ec_stats(new_ec_model)
    
    ec_model = ChU.compressed_model(ec_model)
    new_ec_model = ChU.compressed_model(new_ec_model)
    
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

## ----------------------------------------------------------------------------
# saving
UJL.save_data(ecG.EC_REFERENCE_DATA, refs)


