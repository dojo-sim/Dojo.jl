function get_dzhanibekov(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    springs=0,
    dampers=0, 
    color=RGBA(0.9,0.9,0.9,1),
    limits=false,
    joint_limits=Dict(),
    T=Float64)

    # mechanism
    origin = Origin{Float64}()
    main_body = Capsule(0.1, 1, 1; color, name=:main)
    main_body.inertia = Diagonal([3e-2, 1e-3, 1e-1])
    side_body = Capsule(0.5 * 0.1, 0.35 * 1, 0.5 * 1;
        orientation_offset=RotY(0.5 * π),
        color=color, name=:side)
    bodies = [main_body, side_body]

    joint_float = JointConstraint(Floating(origin, bodies[1]); name=:floating)
    joint_fixed = JointConstraint(Fixed(bodies[1], bodies[2];
        parent_vertex=szeros(3), child_vertex=[-0.25 * 1; 0; 0]), name=:fixed)
    joints = [joint_float, joint_fixed]
    
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
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, joint_float, [0; 0; 1; 0; 0; 0])

    # construction finished
    return mechanism
end

function initialize_dzhanibekov!(mechanism::Mechanism{T,Nn,Ne,Nb};
    linear_velocity=zeros(3),
    angular_velocity=zeros(3)) where {T,Nn,Ne,Nb}

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :floating);
        xmin=[0;0;1;0;0;0; linear_velocity; angular_velocity])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :fixed))
end
