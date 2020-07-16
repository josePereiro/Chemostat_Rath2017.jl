base_intake_info = nothing
function load_base_intake_info()
    !isfile(BASE_INTAKE_INFO_FILE) && return nothing
    global base_intake_info = Dict{String, Dict}()
    df = DataFrame(CSV.read(BASE_INTAKE_INFO_FILE))
    for (id, c, lb) in zip(df[!,1], df[!,2], df[!,3])
        base_intake_info[id] = Dict("c" => c, "lb" => lb)
    end
    return base_intake_info
end
load_base_intake_info()

"""
    returns a copy of the base_intake_info but with the 
    feed medium concentration of a given Rath steady state
"""
function stst_base_intake_info(stst) 
    intake_info = deepcopy(base_intake_info)
    
    # The feed medium of each steady state only vary
    # in the composition of this mets (see Rath2017)
    for rath_met in ["GLC", "GLN", "GAL"]
        model_met = mets_map[rath_met]
        model_exch = exch_met_map[model_met]
        intake_info[model_exch]["c"] = RathData.val("c$rath_met", stst)
    end
    return intake_info
end