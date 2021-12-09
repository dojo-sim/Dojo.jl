# Open visualizer
vis = Visualizer()
open(vis)

################################################################################
# Dzhanibekov effect
################################################################################
Δt0 = 0.01
g0 = 0.0
mech = getmechanism(:dzhanibekov, Δt=Δt0, g=g0)

bodies = collect(mech.bodies)
bodies[1].m
bodies[2].m
bodies[1].J
bodies[2].J

initialize!(mech, :dzhanibekov, x=[0.,0.,0.], v = [0.,0.,0.], q=UnitQuaternion(1.,0.,0.,0.), ω = [5.,0.01,0.])

storage = simulate!(mech, 10.0, record=true)
visualize(mech, storage, vis=vis)