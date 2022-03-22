vis = Visualizer()
open(vis)

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.1, 1.0)
joint0to1 = JointConstraint(Floating(origin, pbody))
bodies = [pbody]
joints = [joint0to1]
contact_p = contact_constraint(pbody, [0.0, 0.0, 1.0], 
    friction_coefficient=0.5, 
    contact_origin=[0.0, 0.0, 0.0], 
    contact_radius=0.1,
    contact_type=:nonlinear)

contacts = [contact_p]
mech = Mechanism(origin, bodies, joints, contacts,
            gravity=-9.81, 
            timestep=0.1)
mech.bodies[1].state.x2 = [0.0, 0.0, 0.5]

storage = simulate!(mech, 1.0, 
    verbose=true, 
    record=true)

visualize(mech, storage, vis=vis)