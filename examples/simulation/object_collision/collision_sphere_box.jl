using Dojo

# Parameters
origin = Origin{Float64}()
pbody = Box(1.0, 1.0, 1.0, 1.0)
cbody = Dojo.Sphere(0.5, 1.0)
joint1 = JointConstraint(Floating(origin, pbody))
joint2 = JointConstraint(Floating(origin, cbody))

bodies = [pbody, cbody]
joints = [joint1, joint2]

side = 1.0
corners = [
            [[ side / 2.0;  side / 2.0; -side / 2.0]]
            [[ side / 2.0; -side / 2.0; -side / 2.0]]
            [[-side / 2.0;  side / 2.0; -side / 2.0]]
            [[-side / 2.0; -side / 2.0; -side / 2.0]]
            [[ side / 2.0;  side / 2.0;  side / 2.0]]
            [[ side / 2.0; -side / 2.0;  side / 2.0]]
            [[-side / 2.0;  side / 2.0;  side / 2.0]]
            [[-side / 2.0; -side / 2.0;  side / 2.0]]
        ]
       
normal = [[0.0, 0.0, 1.0] for i = 1:8]
contact_radius = [0.0 for i = 1:8]
friction_coefficient = 0.5 * ones(8)

box_contacts = contact_constraint(pbody, normal, 
    friction_coefficient=friction_coefficient, 
    contact_origins=corners, 
    contact_radius=contact_radius, 
    contact_type=:nonlinear)

sphere_contact = contact_constraint(cbody, [0.0; 0.0; 1.0], 
    friction_coefficient=0.5, 
    contact_origin=SA[0.0; 0.0; 0.0], 
    contact_radius=0.5)

collision = SphereBoxCollision{Float64,2,3,6}(
    szeros(3),
    1.0,
    1.0,
    2 * 1.0,
    0.5,
)

friction_parameterization = SA{Float64}[
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [
        # corner_contact1, 
        # corner_contact2, 
        # corner_contact3, 
        # corner_contact4, 
        # corner_contact5, 
        # corner_contact6, 
        # corner_contact7, 
        # corner_contact8, 
        box_contacts...,
        sphere_contact,
        ContactConstraint((body_body_contact, cbody.id, pbody.id), name=:body_body)
]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.01)

mech.bodies[1].state.x2 = [0.0, 0.0, 1.5]
q = rand(4) 
q ./= norm(q) 
mech.bodies[1].state.q2 = Quaternion(q...)
mech.bodies[2].state.x2 = [0.3, 0.25, 3.0]

storage = simulate!(mech, 5.0, 
    verbose=true, 
    record=true,
    opts=SolverOptions(verbose=true, rtol=1.0e-5, btol=1.0e-5))

vis = Visualizer()
visualize(mech, storage, 
    vis=vis)
render(vis)
