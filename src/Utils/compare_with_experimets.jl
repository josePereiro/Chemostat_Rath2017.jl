function compare_with_experimets(base_model)
    for stst in RathData.STSTS
        model = deepcopy(base_model)
        println("\nStst: $stst")
        ξ = RathData.val(:ξ, stst)
        println("exp xi: $ξ")
        exp_μ = RathData.val(:μ, stst)
        println("exp growth: $exp_μ")    
        base_intake_info = Human1.HumanGEM.stst_base_intake_info(stst) 
        Chemostat.apply_bound!(model, ξ, base_intake_info; 
            emptyfirst = true, ignore_miss = true    
        )
        fbaout_μ = MetLP.objval(MetLP.fba!(model, Human1.HumanGEM.HUMAN_BIOMASS_IDER))
        println("fba growth: $(fbaout_μ)") 
        fbaout_μ < exp_μ && @warn("fbaout_μ ($fbaout_μ) > exp_μ ($exp_μ)")
    end
end