using Dojo
using FiniteDiff

# Parameters
origin = Origin{Float64}()
pbody = Box(1.0, 1.0, 1.0, 1.0)
cbody = Sphere(0.5, 1.0)
joint = JointConstraint(Fixed(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

collision = SphereSphereCollision{Float64,2,3,6}(
        szeros(3),
        szeros(3),
        # pbody.shape.r, 
        0.5,
        0.5,
        # 0.5,
        )

collision = SphereCapsuleCollision{Float64,2,3,6}(
    szeros(3),
    SA[0.0; 1.0; 0.0],
    SA[0.0; -1.0; 0.0],
    0.5,
    0.5,
    )

collision = SphereBoxCollision{Float64,2,3,6}(
    szeros(3),
    SA[0.5; 0.0; 0.0],
    SA[0.0; 0.5; 0.0],
    SA[0.0; 0.0; 0.5],
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
mech.bodies[2].state.x2 = [0.0, 1.0, 2.0]
mech.bodies[2].state.v15 = [0.0; -1.0; -2.0]

storage = simulate!(mech, 1.0, 
    verbose=true, 
    record=true)

vis = Visualizer()
visualize(mech, storage, 
    vis=vis)
open(vis)

