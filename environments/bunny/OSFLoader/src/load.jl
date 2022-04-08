function load_density_script(; mode::Symbol=:cpu, osf_path=OSF_PATH)
    pushfirst!(pyimport("sys")."path", OSF_PATH)
    if mode == :cpu
        @pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
    else
        @pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
    end
    return nothing
end
