function get_dzhanibekov(;
    timestep=0.01,
    gravity=-9.81,
    color=RGBA(0.9,0.9,0.9,1),
    T=Float64)

    radius = 0.1
    body_length = 1.0
    body_mass = 1.0
    origin = Origin{Float64}()
    main_body = Capsule(radius, body_length, body_mass,
        color=color, name=:main)
    side_body = Capsule(0.5 * radius, 0.35 * body_length, 0.5 * body_mass,
        orientation_offset=RotY(0.5 * π),
        color=color, name=:side)
    links = [main_body, side_body]

    # Joint Constraints
    joint_float = JointConstraint(Floating(origin, links[1]), name=:floating)
    joint_fixed = JointConstraint(Fixed(links[1], links[2];
        parent_vertex=szeros(3),
        child_vertex=[-0.25 * body_length; 0.0; 0.0]), name=:fixed)
    joints = [joint_float, joint_fixed]

    links[1].inertia = Diagonal([3e-2, 1e-3, 1e-1])
    # Mechanism
    return Mechanism(origin, links, joints;
        gravity,
        timestep)
end

function initialize_dzhanibekov!(mechanism::Mechanism{T,Nn,Ne,Nb};
    linear_velocity=zeros(3),
    angular_velocity=zeros(3)) where {T,Nn,Ne,Nb}

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :floating);
        xmin=[0;0;1;0;0;0; linear_velocity; angular_velocity])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :fixed))
end
