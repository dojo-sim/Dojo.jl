function osf_loader_dir()
    @__DIR__
end

using Pkg
Pkg.activate(osf_loader_dir())
using OSFLoader
