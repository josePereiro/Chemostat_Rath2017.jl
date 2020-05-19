
# +
# This script take the data from data/raw, fetched from
# Rath 2017 (https://pure.mpg.de/pubman/item/item_2508673_4)
# and convert it for compatibility with further modeling
# -

# ---
# ### Standard Medium 42 MAX-UB
# ---

# Standard medium original file
# This is the base of the feed medium
stdm_orig = DataFrame(CSV.read(RATH_STDM_ORIG_FILE, delim = "\t"));

# +
# Converting all conc to mM
stdm_conv = DataFrame(stdm_orig);

# Converting all needed to mM
# C(g/L) / MM(g/mol) * 1e3 = mM
idx = findfirst(stdm_conv.id .== "HCO3")
stdm_conv[idx, [:conc, :unit]] .= [(stdm_orig[idx, :conc] / 84) * 1e3, "mM"]
idx = findfirst(stdm_conv.id .== "GAL")
stdm_conv[idx, [:conc, :unit]] .= [(stdm_orig[idx, :conc] / 180) * 1e3, "mM"]
idx = findfirst(stdm_conv.id .== "GLC")
stdm_conv[idx, [:conc, :unit]] .= [(stdm_orig[idx, :conc] / 180) * 1e3, "mM"]

# pH = -log(C(H+)), C(H+) = 10^(-pH) * 1e3 = mM
idx = findfirst(stdm_conv.id .== "H")
stdm_conv[idx, [:conc, :unit]] .= [10.0^(-stdm_orig[idx, :conc]) * 1e3, "mM"]

# C(μM) * 1e-3 = mM
for met in ["GLU", "ALA", "ARG", "ASN", "ASP", "CYS", "GLY", 
            "HIS", "ILE", "LEU", "LAC", "LYS", "MET", "PHE", 
            "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    idx = findfirst(stdm_conv.id .== met)
    stdm_conv[idx, [:conc, :unit]] .= [stdm_orig[idx, :conc] * 1e-3, "mM"]
end
# -

# Saving
if !isfile(RATH_STDM_CONV_FILE)
    CSV.write(RATH_STDM_CONV_FILE, stdm_conv, delim = "\t")
    println("created $(relpath(RATH_STDM_CONV_FILE, PROJ_ROOT))!!!")
end

# ---
# ### Cont cult data
# ---

# ##### Data related names

# measured metabolites
# medium concentration
med_s = ["s$met" for met in msd_mets]

# Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
ρ = 0.25

# #### Data Files

# +
# Original files
cul_data_origs = Dict()
cul_data_convs = Dict()
CUL_DATA_ORIG_FILES = Dict()
CUL_DATA_CONV_FILES = Dict()
for exp in exps
    
    # file names
    filename = "$(RATH_CONT_CUL_DATA_FILE_SUFFIX)_$(exp).tsv"
    
    orig_file_path = joinpath(RATH_RAW_DATA_DIR, filename)
    
    if !isfile(orig_file_path)
        error("$(orig_file_path) not found!!!")
    end
    
    # Original data
    orig_data = DataFrame(CSV.read(orig_file_path, delim = "\t"))
    
    # Converted Data
    conv_data = DataFrame([String, Float64, Float64, String], ["id", "val", "err", "unit"])
    
    #TODO Delete
    data = orig_data


    # Converting cell volume per media volume (CVv) to Cell mass density (gDW/ L)
    # μL/ mL * 1e9     = μm^3/ mL                    (1 μL = 1e9 μm^3)
    # μm^3/ mL * ρ     = pgDW/ mL                    (ρ = 0.25 pgDW/ μm^3) Niklas(2011)
    # pgDW/ mL * 1e-12 = gDW/ mL                     (1 pg = 1 g * 1e-12)
    # gDW/ mL * 1e3    = gDW/ L                      (1 mL = 1e-3 L)
    # Xv = cvrath * 1e9 * ρ * 1e-12 * 1e3
    # Xv = cvrath * ρ (gDW/ L)
    rowidx = findfirst(orig_data.id .== "CVv")
    val = orig_data[rowidx, :val] * ρ
    # Expres error relative to the value
    err = abs(val * orig_data[rowidx, "%err"] / 100)
    unit = "gDW/ L"
    push!(conv_data, ["Xv", val, err, unit])
    
    # This values was already in the desired units
    for ider in [med_s; feed_c; ["D", "μ"]]
        rowidx = findfirst(orig_data.id .== ider)
        
        val = orig_data[rowidx, :val]
        err = abs(val * orig_data[rowidx, "%err"] / 100)
        unit = orig_data[rowidx, :unit]
        push!(conv_data, [ider, val, err, unit])
    end
    
    # Converting fluxes (nmol/ μL hr) to (mmol/ gDW hr)
    # nmol/ μL hr * 1e-6     = mmol/ μL hr          (1 nmol  = 1 mmol * 1e-6)
    # mmol/ μL hr * 1e-9     = mmol/ μm^3 hr        (1 μL = 1e9 μm^3)
    # mmol/ μm^3 hr * (1/ ρ) = mmol/ pgDW hr        (ρ = 0.25 pgDW/ μm^3) Niklas
    # mmol/ pgDW hr * 1e12   = mmol/ gDW hr         (1 pg = 1 g * 1e-12)
    # q = qrath * 1e-6 * 1e-9 * (1/ ρ) * 1e12
    # q = qrath * (1/ ρ) * 1e-3 (mmol/ gDW hr)
    for ider in rate_q
        rowidx = findfirst(orig_data.id .== ider)
        
        val = (orig_data[rowidx, :val] * 1e-3) / ρ
        err = abs(val * orig_data[rowidx, "%err"] / 100)
        unit = "mmol/ gDW hr"
        push!(conv_data, [ider, val, err, unit])
    end
    
    # Converting A1AT flux
    rowidx = findfirst(orig_data.id .== "CD")
    CD = orig_data[rowidx, :val] # (μm)
    CV = (4/3)* π *(CD/2)^3 # (μm³)
    CDW = ρ * CV * 1e-12 # (pgDW) -> (gDW)

    a1at_mw = 54; # kDa (kg/mol)
    rowidx = findfirst(orig_data.id .== "qA1AT")
    qa1at = orig_data[rowidx, :val]/ CDW # (pg/cell d) -> (pg/gDW d)
    qa1at = (qa1at/a1at_mw) * 1e-15 # (pg/gDW d) -> (mol/gDW d)
    qa1at = qa1at * 1e3 # (mol/gDW d) -> (mmol/gDW d)
    qa1at = qa1at/24 # (pg/gDW d) -> (mmol/gDW hr)
    err = abs(qa1at * orig_data[rowidx, "%err"] / 100)
    unit = "mmol/ gDW hr"
    push!(conv_data, ["qA1AT", qa1at, err, unit])
        
    # Saving
    conv_file_path = joinpath(RATH_PROCESSED_DATA_DIR, filename)
    if !isfile(conv_file_path)
        CSV.write(conv_file_path, conv_data, delim = "\t")
        println("created $(relpath(conv_file_path, PROJ_ROOT))!!!")
    end

    cul_data_origs[exp] = orig_data
    cul_data_convs[exp] = conv_data
    CUL_DATA_ORIG_FILES[exp] = orig_file_path
    CUL_DATA_CONV_FILES[exp] = conv_file_path

end

# invitro max fluxes
RATH_MAX_FLUX_ORIG_FILE = joinpath(RATH_RAW_DATA_DIR, "rath2017___max_invitro_fluxs.tsv")

max_flux_data_orig = DataFrame(CSV.read(RATH_MAX_FLUX_ORIG_FILE, delim = "\t"))
max_flux_data_conv = DataFrame(max_flux_data_orig)

# Converting fluxes (nmol/ μL hr) to (mmol/ gDW hr)
# nmol/ μL hr * 1e-6     = mmol/ μL hr          (1 nmol  = 1 mmol * 1e-6)
# mmol/ μL hr * 1e-9     = mmol/ μm^3 hr        (1 μL = 1e9 μm^3)
# mmol/ μm^3 hr * (1/ ρ) = mmol/ pgDW hr        (ρ = 0.25 pgDW/ μm^3) Niklas
# mmol/ pgDW hr * 1e12   = mmol/ gDW hr         (1 pg = 1 g * 1e-12)
# q = qrath * 1e-6 * 1e-9 * (1/ ρ) * 1e12
# q = qrath * (1/ ρ) * 1e-3 (mmol/ gDW hr)
for (enz_i, enz) in enumerate(max_flux_data_orig.id)
    # flux data
    for data_id in ["cul_1", "cul_2", "mean", "mean_err"]
        max_flux_data_conv[enz_i, data_id] = (max_flux_data_orig[enz_i, data_id] * 1e-3) / ρ
    end
end
# unit 
max_flux_data_conv[!, :flux_unit] .= "mmol/ gDW hr";

# Saving
RATH_MAX_FLUX_CONV_FILE = joinpath(RATH_PROCESSED_DATA_DIR, basename(RATH_MAX_FLUX_ORIG_FILE))
if !isfile(RATH_MAX_FLUX_CONV_FILE)
    CSV.write(RATH_MAX_FLUX_CONV_FILE, max_flux_data_conv, delim = "\t")
    println("created $(relpath(RATH_MAX_FLUX_CONV_FILE, PROJ_ROOT))!!!")
end

# -


