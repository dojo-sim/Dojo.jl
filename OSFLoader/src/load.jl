function load_python_script(; mode::Symbol=:cpu, osf_path=OSF_PATH)
    pushfirst!(pyimport("sys")."path", OSF_PATH)
    if mode == :cpu
        @pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
    elseif mode == :gpu
        @pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
    end
    return nothing
end
