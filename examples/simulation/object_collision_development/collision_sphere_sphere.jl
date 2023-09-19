using Dojo
using DojoEnvironments: Z_AXIS
using LinearAlgebra

# Parameters
origin = Origin{Float64}()
pbody = Dojo.Sphere(0.5, 1.0)
cbody = Dojo.Sphere(0.5, 1.0)
joint = JointConstraint(Fixed(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

sphere_contact = ContactConstraint(NonlinearContact(cbody, Z_AXIS, 1.0; contact_radius=0.5))

collision = SphereSphereCollision{Float64,2,3,6}(
        szeros(3), szeros(3), pbody.shape.r, cbody.shape.r)

friction_parameterization = [
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [
    sphere_contact;
    ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)
]

mech = Mechanism(origin, bodies, joints, contacts;
            gravity=-9.81, timestep=0.1)

mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
mech.bodies[2].state.x2 = [randn(2)/10..., 2.0]


storage = simulate!(mech, 3.0; record=true)

visualize(mech, storage)