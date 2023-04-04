function get_cartpole(; 
    timestep=0.1, 
    gravity=[0.0; 0.0; -9.81], 
    spring=0.0, 
    damper=0.0,
    slider_mass=1.0,
    pendulum_mass=1.0,
    pendulum_length=1.0,
    T=Float64)

    # Parameters
    slider_axis = [0.0; 1.0; 0.0]
    pendulum_axis = [1.0; 0.0; 0.0]
    slider_length = 1.0
    radius = 0.075

    # Links
    origin = Origin{Float64}()
    slider = Capsule(1.5 * radius, slider_length, slider_mass, 
        orientation_offset=RotX(0.5 * π), 
        color=RGBA(0.7, 0.7, 0.7, 1.0))
    pendulum = Capsule(radius, pendulum_length, pendulum_mass, 
        color=RGBA(0.7, 0.7, 0.7, 1.0))
    links = [slider, pendulum]

    # Joint Constraints
    joint_origin_slider = JointConstraint(Prismatic(origin, slider, slider_axis;
        spring=spring, damper=damper))
    joint_slider_pendulum = JointConstraint(Revolute(slider, pendulum, pendulum_axis; 
        child_vertex=[0.0; 0.0; -0.5 * pendulum_length],
        spring=spring, damper=damper))
    joints = [joint_origin_slider, joint_slider_pendulum]

    # Mechanism
    mech = Mechanism(origin, links, joints;
        gravity, 
        timestep)

    return mech
end

function initialize_cartpole!(mech::Mechanism{T,Nn,Ne,Nb}; 
    mode=:up, 
    pendulum_length=1.0) where {T,Nn,Ne,Nb}

# origin to slider
set_maximal_configurations!(mech.origin, mech.bodies[1])
set_maximal_velocities!(mech.bodies[1], 
    v=[0.0; 0.0; 0.0],
    ω=zeros(3))

# slider to pendulum
if mode == :down
    set_maximal_configurations!(mech.bodies[1], mech.bodies[2], 
        Δx=[0.0; 0.0; -0.5 * pendulum_length], 
        Δq=RotX(π))
    set_maximal_velocities!(mech.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
elseif mode == :up
    set_maximal_configurations!(mech.bodies[1], mech.bodies[2], 
        Δx=[0.0; 0.0; 0.5 * pendulum_length], 
        Δq=RotX(0))
    set_maximal_velocities!(mech.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
end
end

function mujoco_inertia!(mech)
mech.bodies[1].m = 1.0063
mech.bodies[1].J = Diagonal([0.106974, 0.106974, 0.00636812])

mech.bodies[2].m = 0.4321
mech.bodies[2].J = Diagonal([0.0422274, 0.0422274, 0.0012155])
end