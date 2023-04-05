function get_tennisracket(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    scale=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0, 
    limits=false,
    joint_limits=Dict(),
    T=Float64) 

    origin = Origin{T}(name=:origin)
    bodies = [Box(scale / 25, scale / 2, scale, mass; color, name=:box)]
    joints = [JointConstraint(Floating(origin, bodies[1]); name=:floating_joint)]

    mechanism = Mechanism(origin, bodies, joints; 
        timestep, gravity, input_scaling)

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
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, joints[1], [0, 0, scale/1.9, 0, 0, 0])

    # construction finished
    return mechanism
end

function initialize_tennisracket!(mechanism::Mechanism{T}; 
    x=zeros(3),
    q=one(Quaternion), 
    v=zeros(3),
    ω=zeros(3)) where T

    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [x; rotation_vector(q)])
    set_minimal_velocities!(mechanism, joint, [v; ω])
end
