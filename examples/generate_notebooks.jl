using Literate

function preprocess(str)
    str = replace(str, "# PREAMBLE" => "")
    str = replace(str, "# PKG_SETUP" =>
    """
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    """)
    return str
end

exampledir = @__DIR__
for subdir in ["simulation", "trajectory_optimization", "reinforcement_learning", "gradients", "mechanics"]
    root = joinpath(exampledir, subdir)
    isdir(root) || continue
    @show subdir
    for file in readdir(root)
        name, ext = splitext(file)
        lowercase(ext) == ".jl" || continue
        absfile = joinpath(root, file)
        @show absfile
        Literate.notebook(absfile, root, execute=false, preprocess=preprocess)
    end
end