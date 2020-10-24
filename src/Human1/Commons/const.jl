# A few common useful constants
const ATPM_IDER = "HMR_3964"
const PROT_POOL_EXCHANGE = "prot_pool_exchange"
const PROT_POOL = "prot_pool"
const BIOMASS_IDER = "biomass_human";
const EXCH_SUBSYS_HINT = "Exchange/demand"
const ZEROTH = 1e-8
const MAX_BOUND = 1e3
const UP_FREC = 100
const REV_SUFFIX = "_REV"
const NUM_SUFFIX = "No"
const PROT_PREFFIX = "prot_"
const ARM_PREFFIX = "arm_"
const PMET_PREFFIX = "pmet_"
const DRAW_PREFFIX = "draw_"
const EMPTY_SPOT = "";
const EXTRAS_KEY = "EXTRAS"
const ECMAP_KEY = "ecmap"
const PROTLESS_KEY = "protless"
const EC_MODEL_KEY = "ec_model"
const SRC_MODEL_KEY = "src_model"
const PROT_STOIS_KEY = "prot_stois"

# Cell density ρ = 0.25 pgDW/μm³ from Niklas(2011) https://doi.org/10.1007/s00449-010-0502-y.
# pgDW/μm³ * 1e9 = pgDW/μL
# pgDW/μL * 1e6 = pgDW/L
# pgDW/L * 1e-12 = gDW/L
# atpm_flux = 0.3 mol ATP/ L hr from Fernandez-de-Cossio-Diaz,
# Jorge, and Alexei Vazquez. https://doi.org/10.1038/s41598-018-26724-7.
ATPM_FLUX = 0.3 / (0.25 * 1e9 * 1e6 * 1e-12) * 1e3  # mmol/gWD hr
