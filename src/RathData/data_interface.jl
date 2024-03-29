
# This interface the data stored in data/processed/rath2017__data fetched from 
# Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)


# A single structure with all the data.
const _RATH_BUNDLE = Dict()
function _load_rath_bundle()
    
    files = [_cont_culture_proc_file(exp) for exp in EXPS]
    push!(files, _stand_medium_proc_file())
    for file in files
        !isfile(file) && return _RATH_BUNDLE
    end

    empty!(_RATH_BUNDLE)
    for exp in EXPS
        
        data = get!(_RATH_BUNDLE, exp, Dict())
        
        # std medium
        stdm_conv = CSV.read(_stand_medium_proc_file(), DataFrame);
        for (i, met) in enumerate(stdm_conv[!, :id])
            val = stdm_conv[i, :conc]
            err = 0.0
            unit = stdm_conv[i, :unit]
            # id,                   val,             err,      unit
            data["c$(met)"] = Dict("val" => val, "err" => err, "unit" => unit)
        end
        
        #cont cul data
        cul_data_convs = CSV.read(_cont_culture_proc_file(exp), DataFrame)
        for (i, id) in enumerate(cul_data_convs[!, :id])
            val = cul_data_convs[i, :val]
            err = cul_data_convs[i, :err]
            unit = cul_data_convs[i, :unit]
            # id,       val, err, unit
            data[id] = Dict("val" => val, "err" => err, "unit" => unit)
        end
        
        # ξ
        ξval = data["Xv"]["val"] / data["D"]["val"]
        ξerr = ξval * (data["Xv"]["err"]/ data["Xv"]["val"])
        ξunit = "gDW * hr/ L"
        data["ξ"] = Dict("val" => ξval, "err" => ξerr, "unit" => ξunit)
            
    end
    return _RATH_BUNDLE
end

# interface
# Will define the val, qval, sval, sval, err and unit function to access Rath data
_parse_id(prefix, id) = (string(prefix) == "q" && string(id) == "μ") ? "μ" : string(prefix, id) # handle biomass
function _define_interface()
    for base_fun in [:val, :err, :unit]      
        key = string(base_fun)
        @eval begin
            $base_fun(id, exp) = _RATH_BUNDLE[string(exp)][string(id)][$key]
            function $base_fun(id, exp::AbstractString, deflt)
                try 
                    $base_fun(id, exp); 
                catch err
                    err isa KeyError && return deflt
                    rethrow(err)
                end
            end
            $base_fun(id, exps::Vector, deflt) = [$base_fun(id, exp, deflt) for exp in exps]
            $base_fun(id, exps::Vector) = [$base_fun(id, exp) for exp in exps]
            $base_fun(id) = $base_fun(id, $EXPS) # all the experiments
        end 

        for prefix in ["q", "c", "s"]
            fun = Symbol(string(prefix, base_fun))
            @eval $fun(id) = $base_fun(_parse_id($prefix, id))
            @eval $fun(id, exp) = $base_fun(_parse_id($prefix, id), exp)
            @eval $fun(id, exps, deflt) = $base_fun(_parse_id($prefix, id), exps, deflt)
        end
    end
end