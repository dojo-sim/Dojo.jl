function get_nslider(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    num_bodies=5,
    color=RGBA(1, 0, 0),
    springs=0, 
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Cylinder(0.05, 1, 1; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Prismatic(origin, bodies[1], Z_AXIS))
    joints = JointConstraint{T}[
        jointb1;
        [JointConstraint(Prismatic(bodies[i - 1], bodies[i], Z_AXIS; 
            parent_vertex=-Y_AXIS*0.05, child_vertex=Y_AXIS*0.05)) for i = 2:num_bodies
        ]
    ]

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

    # zero coordinates
    initialize_nslider!(mechanism)
    
    # construction finished
    return mechanism
end

function initialize_nslider!(mechanism::Mechanism; 
    position=mechanism.bodies[1].shape.rh[2], velocity=0)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    set_minimal_coordinates!(mechanism, mechanism.joints[1], [position])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [velocity])

    return
end
