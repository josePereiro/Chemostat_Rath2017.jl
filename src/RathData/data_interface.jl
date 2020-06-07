
# This interface the data stored in data/processed/rath2017__data fetched from 
# Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)


# A single structure with all the data. one ring to rule them all!!!
rath_bundle = nothing
function load_rath_bundle()
    !isfile(RATH_STDM_CONV_FILE) && return nothing
    global rath_bundle = Dict()
    for exp in exps
        
        data = Dict()
        
        # std medium
        stdm_conv = DataFrame(CSV.read(RATH_STDM_CONV_FILE));
        for (i, met) in enumerate(stdm_conv[!, :id])
            val = stdm_conv[i, :conc]
            err = 0.0
            unit = stdm_conv[i, :unit]
            # id,                   val,             err,      unit
            data["c$(met)"] = Dict("val" => val, "err" => err, "unit" => unit)
        end
        
        #cont cul data
        cul_data_convs = DataFrame(CSV.read(RATH_CONT_CUL_DATA_CONV_FILES[exp]))
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
            
        
        rath_bundle[exp] = data
    end
    return rath_bundle
end
load_rath_bundle()


# interface
for fun in ["val", "err", "unit"]      

    eval(Meta.parse("""
        $(fun)(id, exp::AbstractString) = rath_bundle[exp][string(id)]["$(fun)"]"""))

    eval(Meta.parse("""
        $(fun)(id, exps::Vector) = 
                [$(fun)(id, exp) for exp in exps]"""))

    eval(Meta.parse("""
        function $(fun)(id, exp::AbstractString, deflt)
                try
                      return $(fun)(id, exp)
                catch KeyError
                      return deflt
                end
        end """))
    
    eval(Meta.parse("""
        $(fun)(id, exps::Vector, deflt) =
                [$(fun)(id, exp, deflt) for exp in exps]"""))
    
end