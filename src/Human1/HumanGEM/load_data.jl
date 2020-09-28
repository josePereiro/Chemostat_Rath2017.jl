function _load_all_data()
    for (var, file) in [(:mets_map, METS_MAP_FILE), 
                        (:exch_met_map, EXCH_MET_MAP_FILE),
                        (:niklas_biomass, NIKLAS_BIOMASS_FILE),
                        (:readable_met_ids_map, BASE_READABLE_MET_IDS_FILE), 
                        (:ham_medium, HAM_MEDIUM_FILE), 
                        (:base_intake_info, BASE_INTAKE_INFO_FILE)]

        @eval global $var = isfile($file) ? load_data($file) : nothing
    end
end