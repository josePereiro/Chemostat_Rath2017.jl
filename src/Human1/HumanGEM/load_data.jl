function load_data(datafile)
    !isfile(datafile) && return nothing
    dat = Dict()
    df = DataFrame(CSV.read(datafile))
    for (k, v) in zip(df[!,1], df[!,2])
        dat[k] = v
    end
    return dat
end

function load_all_data()
    global mets_map = load_data(METS_MAP_FILE)
    global exch_met_map = load_data(EXCH_MET_MAP_FILE)
    global niklas_biomass = load_data(NIKLAS_BIOMASS_FILE)
    global readable_met_ids_map = load_data(BASE_READABLE_MET_IDS_FILE)
    global ham_medium = load_data(HAM_MEDIUM_FILE)
end
load_all_data()