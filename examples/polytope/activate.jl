function module_dir()
    @__DIR__
end

using Pkg
Pkg.activate(module_dir())
using Polytope
