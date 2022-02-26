function get_tennisracket(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81],
    T=Float64) 

    origin = Origin{T}(name="origin")
    mass = 1.0
    r = 0.1
    h = 1.0
    bodies = [Box(h / 25.0, h / 2, h, mass, color=RGBA(1., 0., 0.), name=:box)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]
    mechanism = Mechanism(origin, bodies, joints, 
        timestep=timestep, 
        gravity=gravity)
    return mechanism
end

function initialize_tennisracket!(mechanism::Mechanism{T}; 
    x=zeros(3),
    q=one(UnitQuaternion), 
    v=zeros(3),
    ω=zeros(3)) where T

    joint = get_joint_constraint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [x; rotation_vector(q)])
    set_minimal_velocities!(mechanism, joint, [v; ω])
end
