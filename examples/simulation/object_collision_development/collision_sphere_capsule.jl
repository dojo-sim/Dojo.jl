using Dojo
using DojoEnvironments: Z_AXIS, Y_AXIS
using LinearAlgebra

# Parameters
origin = Origin{Float64}()
pbody = Capsule(0.5, 2.0, 1.0, orientation_offset=Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0.0, 0.0))
cbody = Dojo.Sphere(0.5, 1.0)
joint1 = JointConstraint(Floating(origin, pbody))
joint2 = JointConstraint(Floating(pbody, cbody))

bodies = [pbody, cbody]
joints = [joint1]

capsule_contact1 = ContactConstraint(NonlinearContact(pbody, Z_AXIS, 1.0; 
                    contact_origin=Y_AXIS, 
                    contact_radius=0.5))

capsule_contact2 = ContactConstraint(NonlinearContact(pbody, Z_AXIS, 1.0;
    contact_origin=-Y_AXIS, 
    contact_radius=0.5))

sphere_contact = ContactConstraint(NonlinearContact(cbody, Z_AXIS, 1.0;
    contact_radius=0.5))

collision = SphereCapsuleCollision{Float64,2,3,6}(
    szeros(3), Y_AXIS, -Y_AXIS, 0.5, 0.5,
)

friction_parameterization = [
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [
        capsule_contact1;
        capsule_contact2; 
        sphere_contact;
        ContactConstraint((body_body_contact, cbody.id, pbody.id), name=:body_body)
]

mech = Mechanism(origin, bodies, joints, contacts;
            gravity=-9.81, timestep=0.05)
mech.bodies[1].state.x2 = [0.25, 0.25, 1.5]
q = LinearAlgebra.normalize(rand(4)) 
mech.bodies[1].state.q2 = Quaternion(q...)
mech.bodies[2].state.x2 = [0.0, 0.5, 4.0]

storage = simulate!(mech, 2.5, record=true;
    opts=SolverOptions(rtol=1.0e-5, btol=1.0e-5))

visualize(mech, storage)