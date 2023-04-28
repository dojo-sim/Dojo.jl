using Dojo
using DojoEnvironments: Z_AXIS
using LinearAlgebra

# Parameters
origin = Origin{Float64}()
box = Box(0.2, 1.0, 1.0, 1.0)
left_finger = Dojo.Sphere(0.1, 0.1; color=RGBA(1,0,0,1))
right_finger = Dojo.Sphere(0.1, 0.1; color=RGBA(1,0,0,1))
joint1 = JointConstraint(PlanarAxis(origin, box,[1;0;0]))
joint2 = JointConstraint(FixedOrientation(origin, left_finger))
joint3 = JointConstraint(FixedOrientation(origin, right_finger))

bodies = [box, left_finger, right_finger]
joints = [joint1, joint2, joint3]

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

collision = SphereBoxCollision{Float64,2,3,6}(
    szeros(3), 1.0, 1.0, 2 * 1.0, 0.1
)

friction_parameterization = [
    1.0  0.0
    0.0  1.0
]
body_body_contact = NonlinearContact{Float64,8}(1.5, friction_parameterization, collision)

contacts = [
        ContactConstraint((body_body_contact, left_finger.id, box.id), name=:body_body1)
        ContactConstraint((body_body_contact, right_finger.id, box.id), name=:body_body2)
]

mech = Mechanism(origin, bodies, joints, contacts;
            gravity=-9.81, timestep=0.01)

function velocity_controller!(mechanism,v_des)
    vy, vz = box.state.v15[2:3]
    ω = box.state.ω15[1]
    
    Δvy = v_des[1] - vy
    Δvz = v_des[2] - vz
    Δω = v_des[3] - ω
    left_y = 1.2     +Δω +0.5Δvy
    left_z = 1.2+0.1             +0.5Δvz
    right_y = -1.2   -Δω +0.5Δvy
    right_z = 0+0.1              +0.5Δvz
    set_input!(joint2, [0;left_y;left_z]*9.81)
    set_input!(joint3, [0;right_y;right_z]*9.81)
end

function position_controller!(mechanism,x_des)
    y, z = box.state.x2[2:3]
    θ = box.state.q2.v1*2 # small angle approximation

    v_des = x_des - [y;z;θ]

    velocity_controller!(mechanism, v_des)
end

function controller!(mechanism,k)
    position_controller!(mechanism,[0.2*(1-cos(k/50));0.75+0.2*sin(k/50);0])
end

mech.bodies[1].state.x2 = [0.0, 0.0, 0.75]
mech.bodies[2].state.x2 = [0.0, -0.601, 0.5]
mech.bodies[3].state.x2 = [0.0, 0.601,1.0]

storage = simulate!(mech, 5.0, controller!; record=true)

visualize(mech, storage)
