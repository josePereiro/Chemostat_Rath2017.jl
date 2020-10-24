load_mets_map() = load_data(METS_MAP_FILE; verbose = false)
load_exch_met_map() = load_data(EXCH_MET_MAP_FILE; verbose = false)
load_niklas_biomass() = load_data(NIKLAS_BIOMASS_FILE; verbose = false)
load_readable_met_ids_map() = load_data(BASE_READABLE_MET_IDS_FILE; verbose = false)
load_ham_medium() = load_data(HAM_MEDIUM_FILE; verbose = false)
load_base_intake_info() = load_data(BASE_INTAKE_INFO_FILE; verbose = false)

function load_rath_met_exch_map()

    exch_met_map = load_exch_met_map()
    mets_map = load_mets_map()
    model_iders_to_plot = map(Rd.iders_to_plot) do rath_ider
        model_ider = rath_ider == Rd.growth_ider ? BIOMASS_IDER :
            exch_met_map[mets_map[rath_ider]]
    end

    iders_to_plot_map = Dict()
    for (rath_ider, model_ider) in 
            zip(Rd.iders_to_plot, model_iders_to_plot)
        iders_to_plot_map[rath_ider] = model_ider
        iders_to_plot_map[model_ider] = rath_ider
    end
    return iders_to_plot_map 
end