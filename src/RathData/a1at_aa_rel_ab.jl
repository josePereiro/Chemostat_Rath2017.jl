function load_a1at_aa_rel_ab()
    _A1AT_AA_REL_AB = Dict()
    datfile = _a1at_aa_rel_abundance_proc_file()
    !isfile(datfile) && error("data file missing ", datfile)
    df = CSV.read(_a1at_aa_rel_abundance_proc_file(), DataFrame)
    for (m1, m2) in zip(df[!,1], df[!,2])
        _A1AT_AA_REL_AB[m1] = m2
    end
    return _A1AT_AA_REL_AB
end
