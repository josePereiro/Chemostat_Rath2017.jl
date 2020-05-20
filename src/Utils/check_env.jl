"""
Throw an error if not in the package project enviroment
"""
function check_env()
    project_toml_file = relpath(joinpath(PROJ_ROOT, "Project.toml"))
    err_msg = """Not in the package enviroment, use Pkg.activate("$(project_toml_file)") " *
                "or run the script with 'julia --project' from the project root " *
                "directory $PROJ_ROOT"""

    cur_projt = Base.current_project()
    isnothing(cur_projt) && error(err_msg)
    abspath(dirname(cur_projt)) != abspath(PROJ_ROOT) && error(err_msg)
    return nothing
end