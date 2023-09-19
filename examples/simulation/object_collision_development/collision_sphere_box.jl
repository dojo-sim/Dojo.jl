using Dojo
using DojoEnvironments: Z_AXIS
using LinearAlgebra

# Parameters
origin = Origin{Float64}()
pbody = Box(1.0, 1.0, 1.0, 1.0)
cbody = Dojo.Sphere(0.5, 1.0)
joint1 = JointConstraint(Floating(origin, pbody))
joint2 = JointConstraint(Floating(origin, cbody))

bodies = [pbody, cbody]
joints = [joint1, joint2]

side = 1.0
contact_origins = [
            [[ side / 2.0;  side / 2.0; -side / 2.0]]
            [[ side / 2.0; -side / 2.0; -side / 2.0]]
            [[-side / 2.0;  side / 2.0; -side / 2.0]]
            [[-side / 2.0; -side / 2.0; -side / 2.0]]
            [[ side / 2.0;  side / 2.0;  side / 2.0]]
            [[ side / 2.0; -side / 2.0;  side / 2.0]]
            [[-side / 2.0;  side / 2.0;  side / 2.0]]
            [[-side / 2.0; -side / 2.0;  side / 2.0]]
        ]
       
normals = fill(Z_AXIS,8)
friction_coefficients = fill(0.5,8)

box_contacts = ContactConstraint(NonlinearContact(pbody, normals, 
    friction_coefficients; contact_origins))

sphere_contact = ContactConstraint(NonlinearContact(cbody, Z_AXIS, 
    0.5; contact_radius=0.5))

collision = SphereBoxCollision{Float64,2,3,6}(
    szeros(3), 1.0, 1.0, 2 * 1.0, 0.5
)

friction_parameterization = [
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [
        box_contacts...,
        sphere_contact,
        ContactConstraint((body_body_contact, cbody.id, pbody.id), name=:body_body)
]

mech = Mechanism(origin, bodies, joints, contacts;
            gravity=-9.81, timestep=0.01)

mech.bodies[1].state.x2 = [0.0, 0.0, 1.5]
q = LinearAlgebra.normalize(rand(4)) 
mech.bodies[1].state.q2 = Quaternion(q...)
mech.bodies[2].state.x2 = [0.3, 0.25, 3.0]

storage = simulate!(mech, 5.0; record=true, 
    opts=SolverOptions(rtol=1.0e-5, btol=1.0e-5))

visualize(mech, storage)
