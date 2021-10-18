function load_humangem_mat_model()
    human_file = rawdir(Human1, [ "Human1_Publication_Data_Scripts/tINIT_GEMs/data" ], "HumanGEM.mat")
    mat_model = MAT.matread(human_file)["ihuman"];
    mat_model = MetNets.to_symbol_dict(mat_model)
    MetNets.reshape_mat_dict(mat_model)
end

function load_humangem_raw_model(model_dict = load_humangem_mat_model())
    return MetNets.MetNet(model_dict; reshape = true)
end

function load_humangem_base_model(;uncompress = false)
    base_model = ldat(HumanGEM,  
        "HumanGEM_base_model", ".jls"; 
        verbose = false
    )
    return uncompress ? MetNets.uncompressed_model(base_model) : base_model
end