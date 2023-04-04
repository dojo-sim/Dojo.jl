function get_tennisracket(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    scale=1
    springs=0.0,
    dampers=0.0, 
    limits=false,
    joint_limits=Dict(),
    T=Float64) 

    origin = Origin{T}(name=:origin)
    bodies = [Box(scale / 25.0, scale / 2.0, scale, mass, color=RGBA(1.0, 0.0, 0.0), name=:box)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]

    mech = Mechanism(origin, bodies, joints; 
        timestep, gravity, input_scaling)

    # springs and dampers
    set_springs!(mech.joints, springs)
    set_dampers!(mech.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mech, joint_limits)

        mech = Mechanism(Origin{T}(), mech.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # construction finished
    return mech
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
