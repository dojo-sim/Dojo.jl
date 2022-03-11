vis = Visualizer()
open(vis)

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.1, 1.0)
cbody = Sphere(0.1, 1.0)
joint0to1 = JointConstraint(Floating(origin, pbody))
joint1to2 = JointConstraint(Floating(pbody, cbody))
bodies = [pbody, cbody]
joints = [joint0to1, joint1to2]
contact_p = contact_constraint(pbody, [0.0, 0.0, 1.0], 
    friction_coefficient=0.5, 
    contact_origin=[0.0, 0.0, 0.0], 
    contact_radius=0.1,
    contact_type=:nonlinear)
contact_c = contact_constraint(cbody, [0.0, 0.0, 1.0], 
    friction_coefficient=0.5, 
    contact_origin=[0.0, 0.0, 0.0], 
    contact_radius=0.1,
    contact_type=:nonlinear)
contacts = [contact_p, contact_c]
mech = Mechanism(origin, bodies, joints, contacts,
            gravity=-9.81, 
            timestep=0.1)
mech.bodies[1].state.x2 = [0.0, 0.25, 0.5]
mech.bodies[2].state.x2 = [0.0, -0.25, 1.0]

storage = simulate!(mech, 1.0, 
    verbose=true, 
    record=true)

visualize(mech, storage, vis=vis)