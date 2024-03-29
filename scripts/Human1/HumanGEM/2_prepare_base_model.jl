using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    import Chemostat
    import Chemostat.MetNets
    import Chemostat.MetLP
    const Ch = Chemostat
    
    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
end

## ---------------------------------------------------------------------
## Description

# This script prepare the original HumanGEM model for modeling. It include some modifications to it (see comment in code) extracted from several data sources. 

## ---------------------------------------------------------------------
# Prepare base base_model
# this will be the base model for all the processing

# Mat Model
@time begin
    mat_model = HG.load_humangem_mat_model()
    base_model = HG.load_humangem_raw_model(mat_model)
    base_model = MetNets.uncompressed_model(base_model)
    M, N = size(base_model)
    println("\nLoaded Mat model")
    println("Base Model: ", (M, N))
end

## ---------------------------------------------------------------------
base_model = ChR.prepare_metnet(HG, base_model; inf_medium = true);

## ---------------------------------------------------------------------
# Saving base_model
HG.save_humangem_base_model(base_model)