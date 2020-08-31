RATH_GROWTH_IDER = "μ"
function _ider_maps()
    global RATH_IDERS_TO_PLOT = [RATH_GROWTH_IDER; RathData.msd_mets]
    global MODEL_IDERS_TO_PLOT = map(RATH_IDERS_TO_PLOT) do rath_ider
        model_ider = rath_ider == RATH_GROWTH_IDER ? OBJ_IDER :
            HumanGEM.exch_met_map[HumanGEM.mets_map[rath_ider]]
    end
    global IDERS_MAP = Dict()
    for (rath_ider, model_ider) in 
            zip(RATH_IDERS_TO_PLOT, MODEL_IDERS_TO_PLOT)
        IDERS_MAP[rath_ider] = model_ider
        IDERS_MAP[model_ider] = rath_ider
    end
end