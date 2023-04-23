function get_npendulum(;
    timestep=0.01,
    input_scaling=timestep,
    gravity=-9.81,
    num_bodies=5,
    mass=1,
    link_length=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0,
    joint_limits=Dict(),
    keep_fixed_joints=true, 
    base_joint_type=:Revolute,
    rest_joint_type=:Revolute,
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(0.05, 0.05, link_length, mass; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Dojo.Prototype(base_joint_type, origin, bodies[1], X_AXIS;
        parent_vertex=(link_length+0.1)*Z_AXIS*num_bodies, child_vertex=Z_AXIS*link_length/2))

    joints = JointConstraint{T}[
        jointb1;
        [
            JointConstraint(Dojo.Prototype(rest_joint_type, bodies[i - 1], bodies[i], X_AXIS;
            parent_vertex=-Z_AXIS*link_length/2, child_vertex=Z_AXIS*link_length/2)) for i = 2:num_bodies
        ]
    ]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_npendulum!(mechanism)

    # construction finished
    return mechanism
end

function initialize_npendulum!(mechanism::Mechanism;
    base_angle=π/4)
    
    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, mechanism.joints[1], [base_angle])

    return
end
