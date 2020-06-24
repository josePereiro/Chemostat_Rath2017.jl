mets_map = nothing
function load_mets_map()
    !isfile(METS_MAP_FILE) && return nothing
    global mets_map = Dict()
    df = DataFrame(CSV.read(METS_MAP_FILE))
    for (m1, m2) in zip(df[!,1], df[!,2])
        mets_map[m1] = m2
    end
    return mets_map
end
load_mets_map()

exch_met_map = nothing
function load_exch_met_map()
    !isfile(EXCH_MET_MAP_FILE) && return nothing
    global exch_met_map = Dict()
    df = DataFrame(CSV.read(EXCH_MET_MAP_FILE))
    for (m1, m2) in zip(df[!,1], df[!,2])
        exch_met_map[m1] = m2
    end
    return exch_met_map
end
load_exch_met_map()

niklas_biomass = nothing
function load_niklas_biomass()
    !isfile(NIKLAS_BIOMASS_FILE) && return nothing
    global niklas_biomass = Dict()
    df = DataFrame(CSV.read(NIKLAS_BIOMASS_FILE))
    for (met, y) in zip(df[!,1], df[!,2])
        niklas_biomass[met] = y
    end
    return niklas_biomass
end
load_niklas_biomass()
