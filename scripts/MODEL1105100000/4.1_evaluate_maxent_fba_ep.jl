import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

# ------------------------------------------------------------------
@time begin

    using Distributed
    using Serialization
    using SparseArrays
    using Dates
    import StatsBase: mean

    # custom packages
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSU = Ch.SimulationUtils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    
    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
    const M = ChR.MODEL1105100000

    import UtilsJL
    const UJL = UtilsJL

    using Plots
end

## ------------------------------------------------------------------
# LOADING INDEX
# INDEX[:DFILE, stst, method]
INDEX = ChU.load_data(joinpath(M.MODEL_PROCESSED_DATA_DIR, "4_index.bson"))
get_dat(method, stst) = 
    deserialize(joinpath(ChR.PROJ_ROOT, INDEX[:DFILE, stst, method]))

function load_model(name, stst; compressed = false)
    MINDEX = UJL.load_data(M.MODEL_INDEX_FILE; verbose = false)
    mfile = MINDEX[stst][name]
    model = deserialize(mfile)
    compressed ? model : ChU.uncompressed_model(model)
end
# ------------------------------------------------------------------
ME_BOUNDED = :ME_BOUNDED

# ------------------------------------------------------------------
fileid = "4.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), 
        M.MODEL_FIGURES_DATA_DIR; params...
    )
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
FLX_IDERS = ["GLC", "LAC", "GLN", "NH4", "GAL", "PYR", "GLU", "ALA", "ASP"]

COLOR_POOL = [:orange, :blue, :red, :black, :violet, 
    :gray, :green, :brown, :magenta]

EXP_COLORS = Dict(
    "B"   => :blue,
    "A"   => :blue,
    "C"   => :blue,
    "F01" => :orange,
    "D"   => :red,
    "E"   => :brown,
)

IDER_COLORS = let
    iders = [FLX_IDERS; "μ"]
    colors = Plots.distinguishable_colors(length(iders))
    Dict(ider => color for (ider, color) in zip(iders, colors))
end

## ------------------------------------------------------------------
# BIOMASS CORRELATIOM
let

end
## ------------------------------------------------------------------
# Total Correlation
let
    method = ME_BOUNDED

    met_map = M.load_mets_map()
    exch_met_map = M.load_exch_met_map()

    ep_p = plot(;title = "EP flx correlation", xlabel = "exp", ylabel = "model")
    fba_p = plot(;title = "FBA flx correlation", xlabel = "exp", ylabel = "model")
    vals = []
    for stst in Rd.ststs
        INDEX[:DFILE, stst, method]
    
        dat = get_dat(method, stst)
        epout = dat[:epout]
        fbaout = dat[:fbaout]
        model = dat[:model]

        for Rd_ider in FLX_IDERS
            try
                model_ider = exch_met_map[met_map[Rd_ider]]

                Rd_av = -Rd.qval(Rd_ider, stst) # inverted
                Rd_err = Rd.qerr(Rd_ider, stst)
                
                ep_av = ChU.av(model, epout, model_ider) 
                ep_err = sqrt.(ChU.va(model, epout, model_ider))
                
                fba_av = ChU.av(model, fbaout, model_ider) 

                sps = (;label = "", color = IDER_COLORS[Rd_ider], 
                    m = 8, alpha = 0.8)
                scatter!(ep_p, [Rd_av], [ep_av]; 
                    xerr = [Rd_err], yerr = [ep_err], sps...
                )

                scatter!(fba_p, [fba_av], [ep_av]; 
                    xerr = [Rd_err], sps...
                )

                push!(vals, Rd_av, ep_av, fba_av)
            catch err; @warn("Fail", Rd_ider, stst, err) end
        end
    end

    sort!(vals)
    lps = (;label = "", ls = :dash, color = :black, lw = 3, alpha = 0.5)
    plot!(ep_p, vals, vals; lps...)
    plot!(fba_p, vals, vals; lps...)
    mysavefig(Plots.Plot[ep_p, fba_p], "flx_correlation"; method)
end  

## -------------------------------------------------------------------
# leyends
# TODO fix this...
let
    for (title, colors) in [
            ("exp", EXP_COLORS), 
            ("iders", IDER_COLORS)
        ]
    p = plot(; framestyle = :none)
        scolors = sort(collect(colors); by = (p) -> string(first(p)))
        for (id, color) in scolors
            scatter!(p, [0], [0];
                thickness_scaling = 1,
                color, ms = 8, label = string(id),
                legendfontsize=10, 
                # size = [300, 900],
                # legend = :left
            )
        end
        mysavefig(p, "$(title)_color_legend")
    end
end