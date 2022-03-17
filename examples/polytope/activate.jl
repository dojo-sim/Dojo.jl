function example_dir()
    @__DIR__
end

using Pkg
Pkg.activate(example_dir())
using Polytope
