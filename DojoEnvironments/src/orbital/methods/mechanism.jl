function get_orbital(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    num_bodies=5,
    springs=0.0, 
    dampers=0.0, 
    limits=false,
    joint_limits=Dict(),
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(0.05, 0.05, 1, 1, color=RGBA(1.0, 0.0, 0.0)) for i = 1:num_bodies]
    
    jointb1 = JointConstraint(Fixed(origin, bodies[1]; 
        child_vertex=Z_AXIS*h/2))
    
    joints = [
        jointb1;
        [
            JointConstraint(Orbital(bodies[i - 1], bodies[i], Z_AXIS; 
            parent_vertex=-Z_AXIS*h/2, child_vertex=Z_AXIS*h/2, 
            spring, damper)) for i = 2:num_bodies
        ]
    ]

    mech = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

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

function initialize_orbital!(mechanism::Mechanism{T}; 
    orientation=[pi / 4.0, pi / 8.0]) where T

    pbody = mechanism.bodies[1]
    joint = mechanism.joints[1]
    vert11 = joint.translational.vertices[2]
    vert12 = -vert11

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody, 
        child_vertex=vert11, 
        Δq=RotX(0.0))

    set_minimal_coordinates!(mechanism, mechanism.joints[2], orientation)

    zero_velocity!(mechanism)
end
