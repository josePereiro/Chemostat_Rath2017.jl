# +
using DataFrames
using Serialization
using Dates
using StatsBase

using Plots
pyplot()
# -

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat

import Chemostat_Rath2017: Chemostat, RathData, HumanGEM
const Ch = Chemostat
const Rd = Chemostat_Rath2017.RathData
const HG = Chemostat_Rath2017.HumanGEM

