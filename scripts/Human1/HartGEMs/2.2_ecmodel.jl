using ProjAssistant
@quickactivate "Chemostat_Rath2017"

@time begin

    using Plots
    import GR
    isinteractive() && GR.inline("png")
    using Statistics
    using ProgressMeter
    using Base.Threads

    using JuMP
    using GLPK
    using Clp
    using Cbc

    import Chemostat
    import Chemostat.MetNets
    import Chemostat.MetLP
    import Chemostat.MetEP
    const Ch = Chemostat

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const H1 = ChR.Human1
    const HG = H1.HumanGEM
    const AG = H1.HartGEMs

end

## ---------------------------------------------------------------------
# Load model
let
    @info("Loading model")
    # TODO: use a brain cotextualized network
    fn = joinpath(@__DIR__, "../../../", "data/raw/Human1", "Human1_Publication_Data_Scripts/ec_GEMs/models/HOP62/ecModel_batch.mat")
    @show isfile(fn)
    ldat(AG, "ec_model", ".jls") do
        _ec_model = MetNets.read_mat(fn)
        MetNets.clampfields!(_ec_model, [:lb, :ub]; abs_max = 999999, zeroth = 1e-16)
        objidx = MetNets.rxnindex(_ec_model, HG.HUMAN_BIOMASS_IDER)
        _ec_model = ChR.prepare_metnet(AG, _ec_model; inf_medium = false);
        _ec_model = MetNets.force_dims(_ec_model)
        MetNets.compressed_model(_ec_model)
    end |> size
end

## ---------------------------------------------------------------------
_load_ecmodel() = MetNets.uncompressed_model(ldat(AG, "ec_model", ".jls"))

## ---------------------------------------------------------------------
# let
#     @time model = _load_ecmodel()
#     biomider = HG.HUMAN_BIOMASS_IDER
#     biomidx = MetNets.rxnindex(model, biomider)
#     S = model.S
#     @show MetNets.nzabs_range(S)
#     @show MetNets.nzabs_range(S[:, biomidx])
#     @time MetNets.clampfields!(model, [:S]; zeroth = 1e-4)
#     @show MetNets.nzabs_range(S)
#     @show MetNets.nzabs_range(S[:, biomidx])
# end

## ---------------------------------------------------------------------
# let
#     model = _load_ecmodel()
#     biomider = HG.HUMAN_BIOMASS_IDER
#     biomidx = MetNets.rxnindex(model, biomider)
#     solver = Cbc.Optimizer
#     lp_model = MetLP.build_lp_model(model; solver)
#     ntries = 2

#     println()
#     println("-"^60)
#     println("-"^60)
#     println("-"^60)
#     ts0 = Float64[]
#     fbaout = nothing
#     ts0 = map(1:ntries) do i
#         println(); println(i, " ", "-"^60)
#         @elapsed fbaout = MetLP.fba(lp_model, biomidx)
#     end
#     sol0 = MetLP.av(fbaout)
    
#     println()
#     println("-"^60)
#     println("-"^60)
#     println("-"^60)
#     lp_model = MetLP.build_lp_model(model; solver)
#     x = MetLP._get_vars(lp_model)
#     MetLP.set_start_value(lp_model, sol0)
#     ts1 = map(1:ntries) do i
#         println(); println(i, " ", "-"^60)
#         @elapsed fbaout = MetLP.fba(lp_model, biomidx)
#     end

#     println("\n"^5)
#     println("-"^60)
#     @show mean(ts0)
#     @show std(ts0)
#     println("-"^60)
#     @show mean(ts1)
#     @show std(ts1)
# end

## ---------------------------------------------------------------------
# find minimum medioum with essential glucose
include("2.0.1_find_glc_essential_medium.jl")
let   
    return # Test
    model = _load_ecmodel()

    exglcidx = MetNets.rxnindex(model, HG.HUMAN_GLC_EX_IDER)
    biomidx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER) 
    costidx = MetNets.rxnindex(model, HG.PROT_POOL_EXCHANGE)

    find_glc_min_glc_essential_medium!(model, biomidx, exglcidx;
        clear_caches = false,
        tries = 20_000,
        var_th = 0.1,
        fba = (model_) -> MetLP.fba(model_, biomidx, costidx)
    )
end

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
let
    
end

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
let
    # prepare model
    
    for lmodel in [
        MetNets.toy_model,
        MetNets.ecoli_core_model,
        # _load_ecmodel,
    ]
        model = lmodel()
        
        @info("Computing")
        # clamp S
        @time Δord0 = MetNets.nzabs_range(model.S)
        @time rank0 = rank(model.S)
        @time cond0 = cond(model.S)
        @time MetNets.clampfields!(model, [:S]; zeroth = 1e-4, abs_max = 1e5)
        @time Δord1 = MetNets.nzabs_range(model.S)
        @time cond1 = cond(model.S)
        @time rank1 = rank(model.S)
        
        @info("Model Info", size = size(model), Δord0, Δord1, rank0, rank1, cond0, cond1)
    end
end

#=
[ Info: Computing toy model
┌ Info: Model Info
│   size = (5, 8)
│   Δord0 = (0.1, 55.0)
│   Δord1 = (0.1, 55.0)
│   rank0 = 5
│   rank1 = 5
│   cond0 = 59.42143629848969
└   cond1 = 59.42143629848969

[ Info: Computing ecoli core
  0.004317 seconds (9 allocations: 102.688 KiB)
┌ Info: Model Info
│   size = (72, 95)
│   Δord0 = (0.0709, 59.81)
│   Δord1 = (0.0709, 59.81)
│   rank0 = 67
│   rank1 = 67
│   cond0 = 5.4836560841253856e17
└   cond1 = 5.4836560841253856e17

[ Info: Computing echuman
┌ Info: Model Info
│   size = (9358, 22695)
│   Δord0 = (2.7778e-11, 77243.0)
│   Δord1 = (0.0001, 77243.0)
│   rank0 = 9127
│   rank1 = 9127
│   cond0 = 3.974762978485945e22
└   cond1 = 2.62455486336173e22

3.487963 seconds (7 allocations: 1.583 GiB, 1.16% gc time)
579.677279 seconds (12 allocations: 1.588 GiB, 0.39% gc time)
574.091478 seconds (12 allocations: 1.588 GiB, 0.70% gc time)
120.969090 seconds (1.91 G allocations: 31.647 GiB, 11.14% gc time)
  3.315092 seconds (7 allocations: 1.583 GiB, 8.46% gc time)
542.414761 seconds (12 allocations: 1.588 GiB, 0.00% gc time)
606.247566 seconds (12 allocations: 1.588 GiB, 0.07% gc time)
=#

## ---------------------------------------------------------------------# plot glucose dependency
## ---------------------------------------------------------------------# plot glucose dependency
## ---------------------------------------------------------------------# plot glucose dependency
include("2.0.3_glc_gradual_limiting_check.jl") 
let
    # prepare model
    # model = _load_ecmodel()
    model = MetNets.ecoli_core_model()
    
    # clamp S
    Δord0 = MetNets.nzabs_range(model.S)
    rank0 = rank(model.S)
    cond0 = cond(model.S)
    MetNets.clampfields!(model, [:S]; zeroth = 1e-4, abs_max = 1e5)
    Δord1 = MetNets.nzabs_range(model.S)
    cond1 = cond(model.S)
    rank1 = rank(model.S)
    
    @info("Model", size = size(model), Δord0, Δord1, rank0, rank1, cond0, cond1)
    return

    exglcidx = MetNets.rxnindex(model, HG.HUMAN_GLC_EX_IDER)
    biomidx = MetNets.rxnindex(model, HG.HUMAN_BIOMASS_IDER) 
    costidx = MetNets.rxnindex(model, HG.PROT_POOL_EXCHANGE)

    for solver_mod in [Clp, Cbc, GLPK]
        solver = solver_mod.Optimizer
        head_hits = ["minimum_medium", "glc_essential"]
        for name in collect(readdir(datdir(AG)))
            heads, params, _ = parse_dfname(name)
            heads != head_hits && continue
            
            println("\n", "-"^60, "\n")
            @show name
            file = joinpath(datdir(AG), name)
            medium_ids = ldat(AG, file)

            p = plot_glc_gradual_limiting_check(
                model, medium_ids, biomidx, exglcidx;
                exglclb0 = -0.5, 
                fba = (model_) -> MetLP.fba(model_, biomidx, costidx),
                fs = [0.0; 10.0 .^ (-3:0.05:0.0)],
                # fs = [0.0, 0.5, 1.0], 
                solver
            )

            sfig(AG, p, heads..., params, (;solver = string(solver_mod)), ".png"; verbose = true)
        end
    end
end