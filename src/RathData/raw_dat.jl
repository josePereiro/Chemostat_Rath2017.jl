
_stand_medium_raw_file() = rawdir(RathData, "rath2017___42_MAX_UB_standard_medium.tsv")

function _exp_medium_raw_file(expid)
    expid = string(expid)
    (expid in ["A", "B", "C"]) && (expid = "ABC")
    rawdir(RathData, "rath2017___feed_medium_exp_$(expid).tsv")
end

function _cont_culture_raw_file(expid)
    expid = string(expid)
    rawdir(RathData, "rath2017___cont_exp_$(expid).tsv")
end

_max_invitro_fluxs_raw_file() = rawdir(RathData, "rath2017___max_invitro_fluxs.tsv")