using ArgParse
using Dates

set = ArgParseSettings()
@add_arg_table! set begin
    "--dry-run"
        help = "Run without consequences, just printing"
        action = :store_true
    "--dir", "-d"
        help = "Specify the source data dir"
        required = true
    "--upwt", "-w"
        help = "The sleep time (sec) between updates"
        default = 60
end
parsed_args = parse_args(set)
dry_run_flag = parsed_args["dry-run"]
src_dir = abspath(parsed_args["dir"])
upwt = parsed_args["upwt"]

!isdir(src_dir) && error("Source dir not found!!!, ", src_dir)
println("Backing up: ", src_dir)

dest_dir = src_dir * "_backup"
while true
    !isdir(dest_dir) && mkpath(dest_dir)
    for file in readdir(src_dir)
        srcfile = joinpath(src_dir, file)
        destfile = joinpath(dest_dir, file)
        (isfile(destfile) && rand() > 0.05) && continue
        !dry_run_flag && cp(srcfile, destfile, force = true)
        println(basename(destfile), " updated!!! size: ", filesize(destfile))
    end
    sleep(upwt)
end