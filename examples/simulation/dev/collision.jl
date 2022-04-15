using Test

# @testset Collision: Sphere-sphere begin 

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.5, 1.0)
cbody = Sphere(0.5, 1.0)
joint = JointConstraint(Fixed(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

collision = SphereSphereCollision{Float64,2,3,6}(
        szeros(3),
        szeros(3),
        pbody.shape.r, 
        cbody.shape.r)
friction_parameterization = SA{Float64}[
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

contacts = [ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.1)

mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]

d = distance(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

@test abs(d - 1.0) < 1.0e-6

cp = contact_point(:parent, mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

@test norm(cp - [0.0; 0.0; 0.5], Inf) < 1.0e-6

cc = contact_point(:child, mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

@test norm(cc - [0.0; 0.0; 1.5], Inf) < 1.0e-6

cn = contact_normal(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

@test norm(cn - [0.0 0.0 -1.0], Inf) < 1.0e-6

ct = contact_tangent(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

@test norm(ct - [0.0 1.0 0.0; -1.0 0.0 0.0], Inf) < 1.0e-6

am = adjacency_matrix(mech.joints, mech.bodies, mech.contacts) 

@test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

storage = simulate!(mech, 2.0, 
    verbose=false, 
    record=true)

vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis)