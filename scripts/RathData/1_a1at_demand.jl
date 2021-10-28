## ----------------------------------------------------------------------------
using ProjAssistant
@quickactivate

@time begin
    import DataFrames: DataFrame
    import CSV

    import Chemostat_Rath2017
    const ChR = Chemostat_Rath2017
    const Rd = ChR.RathData
end

## ----------------------------------------------------------------------------
# The human glycoprotein
# A1AT has a weight of 52 kDa. It consists of a single polypeptide chain with 394
# amino acids, and it has about 14 % carbohydrates by weight.

# from https=>//www.drugbank.ca/polypeptides/P01009
a1at_sec = string(
    "MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQKTDTSHHDQDHPTFNKITPNLAEFAFS",
    "LYRQLAHQSNSTNIFFSPVSIATAFAMLSLGTKADTHDEILEGLNFNLTEIPEAQIHEGF",
    "QELLRTLNQPDSQLQLTTGNGLFLSEGLKLVDKFLEDVKKLYHSEAFTVNFGDTEEAKKQ",
    "INDYVEKGTQGKIVDLVKELDRDTVFALVNYIFFKGKWERPFEVKDTEEEDFHVDQVTTV",
    "KVPMMKRLGMFNIQHCKKLSSWVLLMKYLGNATAIFFLPDEGKLQHLENELTHDIITKFL",
    "ENEDRRSASLHLPKLSITGTYDLKSVLGQLGITKVFSNGADLSGVTEEAPLKLSKAVHKA",
    "VLTIDEKGTEAAGAMFLEAIPMSIPPEVKFNKPFVFLMIEQNTKSPLFMGKVVNPTQK"
)

aa_map = Dict("CYS"=> "C", "ASP"=> "D", "SER"=> "S", "GLN"=> "Q", "LYS"=> "K",
     "ILE"=> "I", "PRO"=> "P", "THR"=> "T", "PHE"=> "F", "ASN"=> "N", 
     "GLY"=> "G", "HIS"=> "H", "LEU"=> "L", "ARG"=> "R", "TRP"=> "W", 
     "ALA"=> "A", "VAL"=>"V", "GLU"=> "E", "TYR"=> "Y", "MET"=> "M")
for (k, v) in aa_map
    aa_map[v] = k
end

# relative abundance
a1at_aa_rel_ab = Dict()
aa_count = length(a1at_sec)
for i in a1at_sec
    aa = aa_map[string(i)]
    if haskey(a1at_aa_rel_ab, aa)
        a1at_aa_rel_ab[aa] += 1 / aa_count
    else
        a1at_aa_rel_ab[aa] = 1 / aa_count
    end
end

## ----------------------------------------------------------------------------s
df = DataFrame(id = collect(keys(a1at_aa_rel_ab)), rel_ab = collect(values(a1at_aa_rel_ab)))
CSV.write(Rd._a1at_aa_rel_abundance_proc_file(), df)
println("created $(relpath(Rd._a1at_aa_rel_abundance_proc_file()))")