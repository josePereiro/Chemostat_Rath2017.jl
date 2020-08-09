# ---
# ## Print functions
# ---

function print_action(wid, pid, state, head, bodyls...; t = Time(now()))
    Core.println("Worker $wid ($pid) $head at $(t) ----------------------------")
    Core.println("\tState: ", state)
    for body in bodyls
        Core.println("\t$body")
    end
    Core.println()
    flush(stdout);
end

print_action(state, head, bodyls...; kwargs...) =
    remotecall_wait(print_action, 1, myid(), getpid(), state, head, bodyls...; kwargs...)

function string_err(err; max_len = 10000)
    s = sprint(showerror, err, catch_backtrace())
    return length(s) > max_len ? s[1:max_len] * "\n..." : s
end


# ---
# ## caching
# ---

const TEMP_CACHE_FILE_PREFIX = "temp_cache"
CACHE_DIR = pwd()

function set_cache_dir(cache_dir)
    !isdir(cache_dir) && error(cache_dir, " not found!!!")
    global CACHE_DIR = cache_dir
end

temp_cache_file(state, cache_dir = CACHE_DIR) = 
    joinpath(cache_dir, "$(TEMP_CACHE_FILE_PREFIX)___$(hash(state)).jls")


function save_cache(data, state, cache_dir = CACHE_DIR)
    tcache_file = temp_cache_file(state, cache_dir) |> relpath
    try
        serialize(tcache_file, data)
    catch err
            print_action(state, "ERROR SAVING CACHE", 
                "cache_file: $tcache_file", 
                "err:        $(string_err(err))")
    end
    print_action(state, "CACHE SAVED", "cache_file: $tcache_file")
end    

function load_cached(state, cache_dir = CACHE_DIR)
    
    tcache_file = temp_cache_file(state, cache_dir) |> relpath
    data = nothing
    if isfile(tcache_file)
        try
            data = deserialize(tcache_file)
        catch err
            print_action(state, "ERROR LOADING CACHE", 
                "cache_file: $tcache_file", 
                "err:        $(string_err(err))")
        end
        print_action(state, "CACHE LOADED", "cache_file: $tcache_file")
    end
    return data
end

function delete_temp_caches(cache_dir = CACHE_DIR; verbose = true)
    tcaches = filter(file -> startswith(file, TEMP_CACHE_FILE_PREFIX), readdir(cache_dir))
    for tc in tcaches
        tc = joinpath(cache_dir, tc)
        rm(tc, force = true)
        verbose && println(relpath(tc), " deleted!!!")
    end
end
