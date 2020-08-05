function _load_data(datafile)
    !isfile(datafile) && return nothing
    dat = Dict()
    df = CSV.read(datafile, DataFrame)
    for (k, v) in zip(df[!,1], df[!,2])
        dat[k] = v
    end
    return dat
end

function _load_all_data()
    global mets_map = _load_data(METS_MAP_FILE)
    global exch_met_map = _load_data(EXCH_MET_MAP_FILE)
    global niklas_biomass = _load_data(NIKLAS_BIOMASS_FILE)
end