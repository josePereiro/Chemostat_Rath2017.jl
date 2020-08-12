# +
using DrWatson
quickactivate(@__DIR__, "Chemostat_Rath2017")

import Chemostat_Rath2017: DATA_KEY, Rep_Human1
const RepH1 = Rep_Human1;

using SparseArrays
import StatsBase: ecdf, mean, std
import Chemostat
import Chemostat.Utils: MetNet


# -

# ---
# ## Loading data

file = RepH1.COMP_FVA_HG_OUTPUT_FILE
file_dat = wload(file);
dat = file_dat[DATA_KEY];
model_syms = [:orig_model, :ec_model];

# ---
# ## Data To Exctract

to_extract = Dict();
for sym in model_syms
    to_extract[sym] = Dict()
end

# ### Diff Comulative Distrubution

function nz_ecdf(model, zeroth = 1e-8, n = 300)
    
    # Comulative
    abs_diff = [abs(ub - lb) for (lb, ub) in zip(model.lb, model.ub)]
    nz_diff = filter(diff -> diff > zeroth, abs_diff)
    cdf = ecdf(nz_diff)
    
    # range
    max_ = min(maximum(model.ub), 1e4)
    min_ = zeroth
    xs = 10.0 .^ range(log10(min_), log10(max_), length = n)
    
    return (xs, cdf.(xs), mean(nz_diff), std(nz_diff))
end

diff_cdf_sym = :diff_cdf
for sym in model_syms
    model = dat[sym];
    xs_, ys_, mean_, std_ = nz_ecdf(model);
    to_extract[sym][diff_cdf_sym] = (xs = xs_, ys = ys_, mean = mean_, std = std_)
end;

# ### All bounds
bound_sym = :bounds
for sym in model_syms
    model = dat[sym];
    to_extract[sym][bound_sym] = (lb = model.lb, ub = model.ub)
end;


# ---
# ## Save

file = RepH1.COMP_FVA_HG_EXTRACTED_DATA_FILE;
tagsave(file, Dict(DATA_KEY => to_extract))
println(relpath(file), " created, size: ", filesize(file), " bytes")


