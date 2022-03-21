vis = Visualizer()
open(vis)

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.5, 1.0)
cbody = Sphere(0.5, 1.0)
joint0to1 = JointConstraint(Fixed(origin, pbody))
# joint0to1 = JointConstraint(Floating(origin, pbody))

# joint1to2 = JointConstraint(Prismatic(pbody, cbody, [0.0; 0.0; 1.0]))
bodies = [pbody, cbody]
joints = [joint0to1]#, joint1to2]

impact_collision = SphereSphereCollision{Float64,0,3,0}(
        szeros(3),
        szeros(3),
        pbody.shape.r, 
        cbody.shape.r)

linear_collision = SphereSphereCollision{Float64,4,3,12}(
        szeros(3),
        szeros(3),
        pbody.shape.r, 
        cbody.shape.r)

nonlinear_collision = SphereSphereCollision{Float64,2,3,6}(
    szeros(3),
    szeros(3),
    pbody.shape.r, 
    cbody.shape.r)

linear_parameterization = SA{Float64}[
        0.0  1.0
        0.0 -1.0
        1.0  0.0
       -1.0  0.0
]

nonlinear_parameterization = SA{Float64}[
        1.0  0.0
        0.0  1.0
]

body_body_contact = ImpactContact{Float64,2}(szeros(Float64, 0, 3), impact_collision)
# body_body_contact = LinearContact{Float64,12}(0.5, linear_parameterization, linear_collision)
# body_body_contact = NonlinearContact{Float64,8}(0.5, nonlinear_parameterization, nonlinear_collision)

contacts = [ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.1)

mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
mech.bodies[2].state.x2 = [0.0, 0.0, 5.0]

# mech = Mechanism(origin, bodies, joints, contacts,
#             gravity=0.0 * -9.81, 
#             timestep=0.01)

# mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
# mech.bodies[2].state.x2 = [2.0, 2.0, 2.0]
# mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
# mech.bodies[2].state.v15 = [-1.0, -1.0, -1.0]

distance(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

contact_point(:parent, mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

contact_point(:child, mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

contact_normal(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

contact_tangent(mech.contacts[1].model.collision, 
    mech.bodies[1].state.x2, mech.bodies[1].state.q2,
    mech.bodies[2].state.x2, mech.bodies[2].state.q2,)

mech.contacts[1].impulses[2]
mech.contacts[1].impulses_dual[2]

@show adjacency_matrix(mech.joints, mech.bodies, mech.contacts) 

storage = simulate!(mech, 2.5, 
    verbose=true, 
    record=true)

visualize(mech, storage, 
    vis=vis)
