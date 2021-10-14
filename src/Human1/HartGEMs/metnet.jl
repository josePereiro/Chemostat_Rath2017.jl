#  models created by contextualizing Human1 using data from Hart et al., “High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities.”

function load_Hart_mat_model(target_model = "GBM")

    hart_file = rawdir(Chemostat_Rath2017, "Human1_Publication_Data_Scripts/tINIT_GEMs/run_tINIT_outputs/Hart2015/tINIT_Hart2015_HumanGEM_outputs.mat")
    dat = MAT.matread(hart_file)["INIT_output"]
    
    target_modeli = findfirst(isequal(target_model), vec(dat["tissues"]))
    return vec(dat["model"])[target_modeli]
    
end

function load_Hart_metnet(target_model = "GBM")
    model_dict = load_Hart_mat_model(target_model)
    return Chemostat.Utils.MetNet(model_dict; reshape = true)
end