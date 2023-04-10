function get_cartpole(; 
    timestep=0.01, 
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
    keep_fixed_joints=false, 
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

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    initialize_cartpole!(mechanism)

    # construction finished
    return mechanism
end

function initialize_cartpole!(mechanism::Mechanism; 
    position=0, orientation=pi/4)

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, mechanism.joints[1], [position])
    set_minimal_coordinates!(mechanism, mechanism.joints[2], [orientation])

    return
end