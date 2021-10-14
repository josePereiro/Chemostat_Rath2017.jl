function load_humangem_mat()
    human_file = rawdir(Chemostat_Rath2017, "Human1_Publication_Data_Scripts/tINIT_GEMs/data/HumanGEM.mat")
    mat_model = MAT.matread(human_file)["ihuman"];
    mat_model = Chemostat.Utils.to_symbol_dict(mat_model)
    Chemostat.Utils.reshape_mat_dict(mat_model)
end

function load_humangem(model_dict = load_humangem_mat())
    return Chemostat.Utils.MetNet(model_dict; reshape = true)
end