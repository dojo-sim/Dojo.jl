function get_cartpole(; 
    timestep=0.1, 
    input_scaling=timestep, 
    gravity=-9.81, 
    slider_mass=1,
    pendulum_mass=1,
    pendulum_length=1,
    radius=0.075,
    color=RGBA(0.7, 0.7, 0.7, 1),
    springs=0, 
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    T=Float64)

    # mechanism
    origin = Origin{Float64}()
    slider = Capsule(1.5 * radius, 1, slider_mass; 
        orientation_offset=RotX(0.5 * π), color)
    pendulum = Capsule(radius, pendulum_length, pendulum_mass; color)
    bodies = [slider, pendulum]
    
    joint_origin_slider = JointConstraint(Prismatic(origin, slider, Y_AXIS))
    joint_slider_pendulum = JointConstraint(Revolute(slider, pendulum, X_AXIS; 
        child_vertex=-0.5*pendulum_length*Z_AXIS))
    joints = [joint_origin_slider, joint_slider_pendulum]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    zero_coordinates!(mechanism)

    # construction finished
    return mechanism
end

function initialize_cartpole!(mech::Mechanism{T,Nn,Ne,Nb}; 
    mode=:up, 
    pendulum_length=1) where {T,Nn,Ne,Nb}

# origin to slider
set_maximal_configurations!(mechanism.origin, mechanism.bodies[1])
set_maximal_velocities!(mechanism.bodies[1], 
    v=[0; 0; 0],
    ω=zeros(3))

# slider to pendulum
if mode == :down
    set_maximal_configurations!(mechanism.bodies[1], mechanism.bodies[2], 
        Δx=[0; 0; -0.5 * pendulum_length], 
        Δq=RotX(π))
    set_maximal_velocities!(mechanism.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
elseif mode == :up
    set_maximal_configurations!(mechanism.bodies[1], mechanism.bodies[2], 
        Δx=[0; 0; 0.5 * pendulum_length], 
        Δq=RotX(0))
    set_maximal_velocities!(mechanism.bodies[2], 
        v=zeros(3), 
        ω=zeros(3))
end
end

function mujoco_inertia!(mechanism)
mechanism.bodies[1].m = 1.0063
mechanism.bodies[1].J = Diagonal([0.106974, 0.106974, 0.00636812])

mechanism.bodies[2].m = 0.4321
mechanism.bodies[2].J = Diagonal([0.0422274, 0.0422274, 0.0012155])
end