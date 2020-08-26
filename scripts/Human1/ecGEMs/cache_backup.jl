src_dir = "/Users/Pereiro/University/Research/Metabolism/MaxEntEP2020/Working_Version/"*
        "Chemostat_Rath2017/data/processed/Human1/ecGEMs/cache";
dest_dir = src_dir * "_backup"
while true
    !isdir(dest_dir) && mkpath(dest_dir)
    for file in readdir(src_dir)
        srcfile = joinpath(src_dir, file)
        destfile = joinpath(dest_dir, file)
        (isfile(destfile) && rand() > 0.05) && continue
        cp(srcfile, destfile, force = true)
        println(basename(destfile), " copied!!! size: ", filesize(destfile))
    end
    sleep(rand(60:120))
end