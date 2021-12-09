# Open visualizer
vis = Visualizer()
open(vis)

################################################################################
# Full Friction Cone vs. Linearized Friction Cone
################################################################################
Δt0 = 0.01
g0 = -9.81
cf0 = 0.30
mech = getmechanism(:box, Δt=Δt0, g=g0, cf=cf0, cone)

bodies = collect(mech.bodies)
bodies[1].m
bodies[1].J

initialize!(mech, :box, x=[0.,0.,1.], v=[4.,2.,1.], q=UnitQuaternion(1.,0.,0.,0.), ω=[0.,0.,0.])
storage = simulate!(mech, 10.0, record=true)
visualize(mech, storage, vis=vis)
