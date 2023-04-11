using Literate

function preprocess(str)
    str = replace(str, "# PKG_SETUP_2" =>
    """
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "../.."))
    Pkg.instantiate()
    """)
    str = replace(str, "# PKG_SETUP" =>
    """
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    """)
    return str
end

exampledir = @__DIR__
for subdir in ["control", "learning", "simulation", "simulation/mechanics", "system_identification"]
    root = joinpath(exampledir, subdir)
    isdir(root) || continue
    @show subdir
    for file in readdir(root)
        name, ext = splitext(file)
        name == "utilities" && continue
        lowercase(ext) == ".jl" || continue
        absfile = joinpath(root, file)
        @show absfile
        Literate.notebook(absfile, root, execute=false, preprocess=preprocess)
    end
end