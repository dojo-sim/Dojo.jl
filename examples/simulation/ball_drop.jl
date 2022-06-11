using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo

timestep=1.0e-2
gravity=-9.81

T=Float64

# Parameters
height = 1.0

# Bodies
origin = Origin{T}()

bodies = [Sphere(0.1, 1.0, name=Symbol("ball"))]

# Joints
joints = [JointConstraint(Floating(origin, bodies[1]))]

# Contacts 
contacts = [contact_constraint(bodies[1], SA[0.0; 0.0; 1.0], 
            friction_coefficient=0.5, 
            contact_origin=SA[0.0; 0.0; 0.0], 
            contact_radius=0.1, 
            contact_type=:nonlinear)]

mech = Mechanism(origin, bodies, joints, contacts,
    gravity=gravity,
    timestep=timestep)

mech.bodies[1].state.x2 = [0.0; 0.0; 0.5]

storage = simulate!(mech, 1.0, 
    record=true, 
    opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6, verbose=true))

vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis)

