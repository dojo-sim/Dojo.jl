using Dojo

# Parameters
origin = Origin{Float64}()
pbody = Capsule(0.5, 2.0, 1.0, axis_offset=Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0.0, 0.0))
cbody = Sphere(0.5, 1.0)
joint = JointConstraint(Floating(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

capsule_contact1 = contact_constraint(pbody, [0.0; 0.0; 1.0], 
                    friction_coefficient=0.5, 
                    contact_origin=SA[0.0; 1.0; 0.0], 
                    contact_radius=0.5)

capsule_contact2 = contact_constraint(pbody, [0.0; 0.0; 1.0], 
    friction_coefficient=0.5, 
    contact_origin=SA[0.0; -1.0; 0.0], 
    contact_radius=0.5)

sphere_contact = contact_constraint(cbody, [0.0; 0.0; 1.0], 
    friction_coefficient=0.5, 
    contact_origin=SA[0.0; 0.0; 0.0], 
    contact_radius=0.5)

collision = SphereCapsuleCollision{Float64,2,3,6}(
    szeros(3),
    SA[0.0; 1.0; 0.0],
    SA[0.0; -1.0; 0.0],
    0.5,
    0.5,
)

friction_parameterization = SA{Float64}[
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [capsule_contact1, capsule_contact2, sphere_contact, ContactConstraint((body_body_contact, cbody.id, pbody.id), name=:body_body)]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.05)
mech.bodies[1].state.x2 = [0.25, 0.25, 1.0]
q = rand(4)
q ./= norm(q)
mech.bodies[1].state.q2 = Quaternion(q..., true)
mech.bodies[2].state.x2 = [0.0, 0.5, 3.0]

storage = simulate!(mech, 5.0, 
    verbose=true, 
    record=true)

vis = Visualizer()
visualize(mech, storage, 
    vis=vis)
open(vis)