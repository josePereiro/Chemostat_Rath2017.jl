a1at_aa_rel_ab = nothing
function load_a1at_aa_rel_ab()
    !isfile(A1AT_AA_REL_ABUNDANCE_FILE) && return nothing
    global a1at_aa_rel_ab = Dict()
    df = DataFrame(CSV.read(A1AT_AA_REL_ABUNDANCE_FILE))
    for (m1, m2) in zip(df[!,1], df[!,2])
        a1at_aa_rel_ab[m1] = m2
    end
    return a1at_aa_rel_ab
end
load_a1at_aa_rel_ab()