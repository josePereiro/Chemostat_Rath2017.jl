function prepare_metnet(proj, base_model; inf_medium = false)

    base_model = deepcopy(base_model)
    
    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    ## Description
    
    # This script prepare the original HumanGEM model for modeling. It include some modifications to it (see comment in code) extracted from several data sources. 
    
    ## ---------------------------------------------------------------------
    make_supp_files = (proj == HG)
    make_supp_files && @info("Creating support files")
    
    ## ---------------------------------------------------------------------
    M, N = size(base_model)
    println("Base Model: ", (M, N))

    ## ---------------------------------------------------------------------
    # mets_readable_ids
    # we build a more human readable met ids
    if make_supp_files
        println("\nMet readable ids")
        mat_model = HG.load_humangem_mat_model()
        met_readable_ids = Dict()
        for i in 1:M
            readable_id = mat_model[:metNames][i] * "[" * mat_model[:comps][mat_model[:metComps][i]] * "]";
            met_readable_ids[base_model.mets[i]] = readable_id
            met_readable_ids[readable_id] = base_model.mets[i]
        end

        # Saving
        sdat(HG, met_readable_ids, 
            "met_readable_ids", ".jls"; 
            verbose = true
        )
    else
        met_readable_ids = HG.load_met_readable_ids()
    end

    ## ---------------------------------------------------------------------
    # Exchanges
    # bkwd and fwd splitted reactions are troublemakers for EP, but they are necessary to model enzymatic costs.
    # So, we leave as least as possible. We unified the exchanges (make them a unique rxn),
    # and let them in a semi-open state (intake bloked, outtake open)

    # I will delete all the Boundary (x) (see comps in the matmodel) metabilites, 
    # leaving only the Extracellular (s) metabolites in the exchange reactions. 
    # Why? they are not required
    println("\nDeleting Boundary (x) metabolites ")
    to_del = [met for met in base_model.mets if endswith(met, "x")];
    println("todel: ", length(to_del))
    println("before: ", size(base_model))
    MetNets.del_met!(base_model, to_del);
    MetNets.del_void_rxns_mets!(base_model)
    base_model = MetNets.compacted_model(base_model)
    println("after: ", size(base_model))

    # redefining the Chemostat exchange criteria
    exch_subsys_hint = "Exchange/demand"
    function is_exchange(model::MetNets.MetNet, ider::MetNets.IDER_TYPE)
        idx = MetNets.rxnindex(model, ider)
        subsys = model.subSystems[idx]
        return occursin(exch_subsys_hint, string(subsys))
    end

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # Exchanges
    exchs = []
    for exch in filter((ider) -> is_exchange(base_model, ider), base_model.rxns)
        exch_i = MetNets.rxnindex(base_model, exch)
        
        # First close it, later if it is what I want, open the outtake
        MetNets.lb!(base_model, exch_i, 0.0)
        MetNets.ub!(base_model, exch_i, 0.0)
        
        mets = MetNets.rxn_mets(base_model, exch_i)
        reacts = MetNets.rxn_reacts(base_model, exch_i)
        
        # I want only the forward monomoleculars
        !(length(mets) == length(reacts) == 1) && continue 
        # I want only the exchanges 's'
        isempty(reacts) && !endswith(base_model.mets[first(reacts)], "s") && continue 
        
        # Because, this reactions are forward unbalanced (A <-> nothing)
        # positibe (+) bounds limit the outtake of the cell and
        # negative (-) bounds limit the intake.
        # Because in the Chemostat the intakes is 
        # controlled by the medium, we'll close all intakes to handle with them later
        # We'll open all outtakes
        MetNets.lb!(base_model, exch_i, 0.0)
        MetNets.ub!(base_model, exch_i, HG.ABS_MAX_BOUND)
        
        push!(exchs, exch_i)
        
    end
    println("\nExchanges: ", exchs |> length)

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # Exch Met map
    # A fast way to get the exch reaction from the metabolite and viceversa
    if make_supp_files
        exch_met_map = Dict()
        for rxn in exchs
            mets = MetNets.rxn_mets(base_model, rxn)
            if length(mets) == 1 # only the monomoleculars
                met = base_model.mets[mets[1]]
                rxn = base_model.rxns[rxn]
                
                exch_met_map[met] = rxn
                exch_met_map[rxn] = met
            end
        end
        println("\nExchanges met map: ", exch_met_map |> length)

        # Saving
        sdat(HG, exch_met_map, 
            "exch_met_map", ".jls"; 
            verbose = true
        )
    else
        exch_met_map = HG.load_exch_met_map()
    end

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # Base intake info
    # The Base model will have a medium (expressed as open intake fluxes) that resamble the cultivation
    # at xi = 1 using the set up in Rath 2017 exp A. So the lb of the intakes will be directly the (negative)
    # concentration in 42_MAX_UB standard medium (see Cossio's paper). Also, a few intakes, not justified in
    # the standard medium will be add based in the intakes of the original model FBA analysis.
    base_intake_info = HG.load_base_intake_info(; inf_medium)

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # Apply medium
    # The steady state assumption in the context of the Chemostat culture impose a constraint over the intakes dependent of xi and c
    ξ = 1
    println("\nApplaying chemostat bound, xi: ", ξ)
    Chemostat.apply_bound!(base_model, ξ, base_intake_info;
        emptyfirst = true, ignore_miss = true
    );

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # Niklas Biomasss

    # I will modified the biomass equation of MODEL1105100000 model with data
    # derived from Niklas (2013): https://doi.org/10.1016/j.ymben.2013.01.002. Table1. (see README)
    # I do not touch the energetic part of the equation, atp + h20 -> adp + h2 + pi
    println("\nApplying Niklas Biomass")

    biomass_idx = MetNets.rxnindex(base_model, HG.HUMAN_BIOMASS_IDER)
    base_model.S[:, biomass_idx] .= zeros(size(base_model, 1))
    for (met, y) in HG.load_Niklas_biomass()
        MetNets.S!(base_model, met, biomass_idx, y)
    end

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # ATPM demand

    # +
    # Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
    # pgDW/μm³ * 1e9 = pgDW/μL
    # pgDW/μL * 1e6 = pgDW/L
    # pgDW/L * 1e-12 = gDW/L
    # atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
    # Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
    atpm_flux = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr

    # This reaction is foward defined with respect to atp
    # atp + more_reacts <-> adp + more_products
    # so we limit the lower bounds as the minimum atp demand 
    # that the cell metabolism must fullfill
    if HG.HUMAN_ATPM_IDER in base_model.rxns
        MetNets.lb!(base_model, HG.HUMAN_ATPM_IDER, atpm_flux)
        
        println("\nATP demand")
        println(MetNets.rxn_str(base_model, HG.HUMAN_ATPM_IDER), " ", MetNets.bounds(base_model, HG.HUMAN_ATPM_IDER))
    end

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # dimentions
    base_model = MetNets.force_dims(base_model)

    ## ---------------------------------------------------------------------
    println("FBOut: ", MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER).obj_val)

    ## ---------------------------------------------------------------------
    # summary
    println("\nBase Model summary")
    MetNets.summary(base_model)

    ## ---------------------------------------------------------------------
    #  FBA Test
    fbaout = MetLP.fba!(base_model, HG.HUMAN_BIOMASS_IDER)
    println("\nFBAout summary")
    MetNets.summary(base_model, fbaout)

    ## ---------------------------------------------------------------------
    println("\nComparing with experiments")
    compare_with_experimets(base_model)

    return base_model
end