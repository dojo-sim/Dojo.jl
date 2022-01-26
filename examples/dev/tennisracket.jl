# Load packages
using MeshCat
using LinearAlgebra

# Open visualizer
vis = Visualizer()
open(vis)

mech = getmechanism(:tennisracket, timestep = 0.01, g = -0.00);
body1 = collect(mech.bodies)[1]
body1.J = Matrix(Diagonal([1,2,3.]))
initialize!(mech, :rectangle, x = [0,0,1.0], q = UnitQuaternion(RotX(0.00*π)), ω = [0,5.0,0.01])
@elapsed storage = simulate!(mech, 8.0, record = true, verbose = false, opts=SolverOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)
