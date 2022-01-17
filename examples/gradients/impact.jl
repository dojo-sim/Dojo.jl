# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


mech = getmechanism(:box2d, Δt=0.01, g=-9.81, contact=true, conetype=:impact);
initialize!(mech, :box2d, x=[0,1.], v=[0,0.], θ=0.0, ω=0.0)

function controller!(mechanism, k)
    setControl!(mech, [0.1,0,0])
    return
end

u0 = [0.0, 0.0, 0.0]
storage = simulate!(mech, 3.00, controller!, record = true, verbose = true,
    opts=InteriorPointOptions(rtol=3e-4, btol=3e-4))
visualize(mech, storage, vis = vis)

z0 = getMaxState(mech)
nu = controldim(mech)
u0 = zeros(nu)
opts_grad = InteriorPointOptions(rtol=1e-8, btol=1e-8)
∇x, ∇u = getMinGradients!(mech, z0, u0, opts=opts_grad)
∇u[:, 1:3]

# apparently, because the impact force completely compensates the gravity force,
# when we compute the gradient of the vertical displacement of the box wrt vertical
# force it is equal to Δt. I.e. we think that applying a veritcal force will
# automatically result in a vertical displacement. In reality, only forces above
# the gravity force will generate vertical displacement.
