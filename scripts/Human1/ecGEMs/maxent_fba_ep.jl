using DrWatson

import MAT
using SparseArrays
using Test
using Distributions

import Chemostat
import Chemostat.Utils: MetNet
const Ch = Chemostat
import Chemostat_Rath2017: DATA_KEY, Human1
const RepH1 = Human1.Rep_Human1;
const ecG = Human1.ecGEMs
const tIG = Human1.tINIT_GEMs;
const HG = Human1.HumanGEM;

## Loading base models