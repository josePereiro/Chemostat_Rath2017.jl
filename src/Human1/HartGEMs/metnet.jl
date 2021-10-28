#  models created by contextualizing Human1 using data from Hart et al., “High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities.”

const HART_GEM_PREFIX = "HartGEM"

function load_Hart_raw_dat()
    hart_file = rawdir(Human1, [ "Human1_Publication_Data_Scripts/tINIT_GEMs/run_tINIT_outputs/Hart2015" ], "tINIT_Hart2015_HumanGEM_outputs.mat")
    return MAT.matread(hart_file)["INIT_output"]
end

function load_Hart_mat_model(;tissue = "GBM")
    dat = load_Hart_raw_dat()
    target_modeli = findfirst(isequal(tissue), vec(dat["tissues"]))
    return vec(dat["model"])[target_modeli]
end

function load_Hart_raw_model(; uncompress = false, tissue = "GBM")
    model_dict = load_Hart_mat_model(;tissue)
    model = MetNets.MetNet(model_dict; reshape = true)
    uncompress ? model : MetNets.compressed_model(model) 
end

save_model(model::MetNets.MetNet; verbose = true, model_params...) = 
    sdat(HartGEMs, MetNets.compressed_model(model), 
        HART_GEM_PREFIX, model_params, ".jls";  verbose
    )

function load_model(;uncompress = false, verbose = false,
        model_params...
    )

    modelid = get(model_params, :modelid, "") |> string
    tissue = get(model_params, :tissue, "") |> string
    stst = get(model_params, :stst, "") |> string

    if modelid == "mat"
        return load_Hart_mat_model(;tissue)
    elseif modelid == "raw"
        return load_Hart_raw_model(;tissue, uncompress)
    elseif modelid == "base" || modelid == "fva_base" || modelid == "scaled" || modelid == "fva_scaled"
        model = ldat(HartGEMs, HART_GEM_PREFIX, (;modelid, tissue), ".jls"; verbose)
        return uncompress ? MetNets.uncompressed_model(model) : model
    elseif modelid == "stst_base" || modelid == "stst_scaled"
        model = ldat(HartGEMs, HART_GEM_PREFIX, (;modelid, tissue, stst), ".jls"; verbose)
        return uncompress ? MetNets.uncompressed_model(model) : model
    else
        error("Unknown '$(modelid)' modelid")
    end
end

function check_modelfile(;params...) 
    exist = isfile(datdir(HartGEMs, HART_GEM_PREFIX, params, ".jls"))
    exist && @info("Model cached ", params...)
    return exist
end