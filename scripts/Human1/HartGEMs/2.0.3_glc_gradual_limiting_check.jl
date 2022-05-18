function plot_glc_gradual_limiting_check(model, medium_ids, objidx, exglcidx;
        exglclb0 = -0.5, 
        fs = [0.0; 10.0 .^ (-3:0.01:0.0)], 
        solver = Clp.Optimizer, 
        fba::Function = (model_) -> MetLP.fba(model_, objidx)
    )

    println("\n", "-"^60, "\n")
    
    # medium ids
    println()
    @info("Minimum medium")
    @show medium_ids
    
    # model exchanges
    println()
    @info("Finding open exchanges")
    exch_idxs = filter(MetNets.exchanges(model)) do exchi
        MetNets.lb(model, exchi) < 0
    end
    exch_ids = model.rxns[exch_idxs]
    @show length(exch_idxs)
    
    # lp_model
    println()
    @info("Building lp_model")
    @time lp_model = MetLP.build_lp_model(model, solver)

    # MetNets.summary.([model], medium_ids)
    # println("\n", "-"^60)
    # @info("Diff")
    # MetNets.summary.([model], setdiff(Set(exch_ids), Set(medium_ids)))
    medium_idxs = MetNets.rxnindex.([model], medium_ids)
    for exchi in exch_idxs
        MetNets.lb!(lp_model, exchi, -0.0)
    end
    for exchid in medium_idxs
        MetNets.lb!(lp_model, exchid, -1000.0)
    end
    exglclb = MetNets.lb(lp_model, exglcidx)
    @show exglclb

    # Set silent
    # JuMP.MOI.set(lp_model, JuMP.MOI.RawParameter("print_level"), 1)

    # full open biomass
    @time begin
        println()
        fbaout = fba(lp_model)
        biom_full = MetNets.av(fbaout, objidx)
        @show biom_full
        println()
    end
    
    bioms_ = []
    fs_ = []
    for f in sort!(collect(fs); rev = true)
        thid = threadid()
        
        glclb = exglclb0 * f
        
        # update lb
        @time begin
            MetLP.lb!(lp_model, exglcidx, glclb)
            fbaout = fba(lp_model)
            biom_glcfree = MetNets.av(fbaout, objidx)

            push!(bioms_, biom_glcfree)
            push!(fs_, f)
            @info("At", thid, f, glclb, biom_glcfree)
        end
    end

    # plot
    p = plot(fs, fill(biom_full, length(fs)); label = "", ls = :dash, c = :red, lw = 3)
    scatter!(p, fs_, bioms_; xlabel = "factor", ylabel = "biom", label = "", c = :blue, ms = 5)
    return p
end