# Open visualizer
vis = Visualizer()
open(vis)

################################################################################
# Dzhanibekov effect
################################################################################
vis = Visualizer()
open(vis)

Δt0 = 0.01
g0 = 0.0
mech = getmechanism(:dzhanibekov, Δt=Δt0, g=g0);

bodies = collect(mech.bodies)
bodies[1].m
bodies[2].m = 0.1
bodies[1].J = Diagonal([1.0, 1.1, 1.0])
bodies[2].J

initialize!(mech, :dzhanibekov, x=[0.,0.,0.], v = [0.,0.,0.], q=UnitQuaternion(1.,0.,0.,0.), ω = [10.0, 0.01,0.0])


storage = simulate!(mech, 5.0, record=true, verbose=false)
visualize(mech, storage, vis=vis)