using Dojo
using FiniteDiff
# @testset Collision: Sphere-sphere begin 

# Parameters
origin = Origin{Float64}()
pbody = Capsule(0.5, 2.0, 1.0, axis_offset=Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0.0, 0.0))
cbody = Sphere(0.5, 1.0)
joint = JointConstraint(Fixed(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

# collision = SphereSphereCollision{Float64,2,3,6}(
#         szeros(3),
#         szeros(3),
#         # pbody.shape.r, 
#         0.5,
#         cbody.shape.r,
#         # 0.5,
#         )

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

contacts = [ContactConstraint((body_body_contact, cbody.id, pbody.id), name=:body_body)]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.05)

mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
mech.bodies[2].state.x2 = [0.2, 0.5, 2.0]

# d = distance(mech.contacts[1].model.collision, 
#     mech.bodies[2].state.x2, mech.bodies[2].state.q2,
#     mech.bodies[1].state.x2, mech.bodies[1].state.q2,)

# @test abs(d - 1.0) < 1.0e-6

# cp = contact_point(:parent, mech.contacts[1].model.collision, 
#     mech.bodies[2].state.x2, mech.bodies[2].state.q2,
#     mech.bodies[1].state.x2, mech.bodies[1].state.q2,)

# @test norm(cp - [0.0; 0.0; 0.5], Inf) < 1.0e-6

# cc = contact_point(:child, mech.contacts[1].model.collision, 
#     mech.bodies[2].state.x2, mech.bodies[2].state.q2,
#     mech.bodies[1].state.x2, mech.bodies[1].state.q2,
#     )

# @test norm(cc - [0.0; 0.0; 1.5], Inf) < 1.0e-6

# cn = contact_normal(mech.contacts[1].model.collision, 
#     mech.bodies[2].state.x2, mech.bodies[2].state.q2,
#     mech.bodies[1].state.x2, mech.bodies[1].state.q2,
#     )

# @test norm(cn - [0.0 0.0 -1.0], Inf) < 1.0e-6

# ct = contact_tangent(mech.contacts[1].model.collision, 
#     mech.bodies[2].state.x2, mech.bodies[2].state.q2,
#     mech.bodies[1].state.x2, mech.bodies[1].state.q2,
#     )

# @test norm(ct - [0.0 1.0 0.0; -1.0 0.0 0.0], Inf) < 1.0e-6

# am = adjacency_matrix(mech.joints, mech.bodies, mech.contacts) 

# @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

storage = simulate!(mech, 1.0, 
    verbose=true, 
    record=true)

vis = Visualizer()
visualize(mech, storage, 
    vis=vis)
open(vis)