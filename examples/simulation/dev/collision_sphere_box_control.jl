using Dojo

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.5, 10.0)
cbody = Box(1.0, 1.0, 2.0, 1.0)
joint1 = JointConstraint(Floating(origin, pbody))
joint2 = JointConstraint(Floating(origin, cbody))

bodies = [pbody, cbody]
joints = [joint1, joint2]

side = 1.0
corners = [
            [[ side / 2.0;  side / 2.0; -side]]
            [[ side / 2.0; -side / 2.0; -side]]
            [[-side / 2.0;  side / 2.0; -side]]
            [[-side / 2.0; -side / 2.0; -side]]
            [[ side / 2.0;  side / 2.0;  side]]
            [[ side / 2.0; -side / 2.0;  side]]
            [[-side / 2.0;  side / 2.0;  side]]
            [[-side / 2.0; -side / 2.0;  side]]
        ]
       
normal = [[0.0, 0.0, 1.0] for i = 1:8]
contact_radius = [0.0 for i = 1:8]
friction_coefficient = 0.5 * ones(8)

sphere_contact = contact_constraint(pbody, [0.0; 0.0; 1.0], 
    friction_coefficient=0.5, 
    contact_origin=SA[0.0; 0.0; 0.0], 
    contact_radius=0.5)

box_contacts = contact_constraint(cbody, normal, 
    friction_coefficient=friction_coefficient, 
    contact_origins=corners, 
    contact_radius=contact_radius, 
    contact_type=:nonlinear)

collision = SphereBoxCollision{Float64,2,3,6}(
    szeros(3),
    side,
    side,
    2 * side,
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
        ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)
]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=0.0 * -9.81, 
            timestep=0.01)

mech.bodies[1].state.x2 = [2.0, 0.1, 2.0]
mech.bodies[1].state.v15 = [-1.0, 0.1, 0.0]
            
mech.bodies[2].state.x2 = [0.0, 0.0, side]
# q = rand(4) 
# q ./= norm(q) 
mech.bodies[2].state.q2 = RotZ(0.25 * π)

function ctrl!(m, k)
    u = [
          0.0; 0.0; 0.0 * -mech.bodies[1].mass * 9.81 * mech.timestep; 0.0; 0.0; 0.0; 0.0; 0.0; -1.0 * mech.bodies[2].mass * 9.81 * mech.timestep; 0.0; 0.0; 0.0;            
    ]
    set_input!(m, u)
    return nothing
end

storage = simulate!(mech, 5.0, 
    ctrl!,
    # verbose=true, 
    record=true,
    opts=SolverOptions(verbose=true, rtol=1.0e-5, btol=1.0e-5));

mech.bodies[1].name
mech.bodies[2].name
vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis)
mech.bodies[2].state.vsol
x31 = next_position(mech.bodies[1].state.x2, mech.bodies[1].state.vsol[1], mech.timestep)
x32 = next_position(mech.bodies[2].state.x2, mech.bodies[2].state.vsol[2], mech.timestep)
q31 = next_orientation(mech.bodies[1].state.q2, mech.bodies[1].state.ϕsol[1], mech.timestep)
q32 = next_orientation(mech.bodies[2].state.q2, mech.bodies[2].state.ϕsol[2], mech.timestep)

distance(collision, x31, q31, x32, q32)
∂distance∂x(:parent, collision, x31, q31, x32, q32)
∂distance∂x(:child, collision, x31, q31, x32, q32)

pp = contact_point(:parent, collision, x31, q31, x32, q32) 
pc = contact_point(:child, collision, x31, q31, x32, q32) 

cop = contact_point_origin(x31, q31, collision.origin_sphere)
coc = contact_point_box(cop, x32, q32, collision.width, collision.depth, collision.height)
