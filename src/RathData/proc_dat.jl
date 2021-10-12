_stand_medium_proc_file() = procdir(RathData, "rath2017___42_MAX_UB_standard_medium.tsv")

_a1at_aa_rel_abundance_proc_file() = procdir(RathData, "a1at_aa_rel_abundance.csv")

function _cont_culture_proc_file(expid)
    expid = string(expid)
    procdir(RathData, "rath2017___cont_exp_$(expid).tsv")
end

_max_invitro_fluxs_proc_file() = procdir(RathData, "rath2017___max_invitro_fluxs.tsv")