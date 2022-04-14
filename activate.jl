function module_dir()
    @__DIR__
end

using Pkg
Pkg.activate(module_dir())
Pkg.instantiate()
using Dojo
