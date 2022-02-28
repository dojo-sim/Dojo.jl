function get_tippetop(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.8, 
    contact=true, 
    contact_type=:nonlinear,
    T=Float64)

    origin = Origin{T}(name=:origin)
    radius = 0.5
    mass = 1.0
    α = 0.2
    bodies = [Sphere(radius, mass, name=:sphere1, color=cyan), Sphere(radius*α, mass*α, name=:sphere2, color=cyan)]
    joints = [JointConstraint(Floating(origin, bodies[1]), 
                name=:floating_joint),
              JointConstraint(Fixed(bodies[1], bodies[2], 
                parent_vertex=[0,0,radius], 
                child_vertex=zeros(3)), name = :fixed_joint),]
    mechanism = Mechanism(origin, bodies, joints,
        timestep=timestep, 
        gravity=gravity)

    # modify inertias
    mechanism.bodies[1].inertia = Diagonal([1.9, 2.1, 2.0])

    if contact
        contact = [0.0, 0.0, 0.0]
        normal = [0.0, 0.0, 1.0]
        contacts = [
            contact_constraint(get_body(mechanism, :sphere1), normal, 
                friction_coefficient=friction_coefficient, 
                contact_point=contact, offset=[0.0, 0.0, radius], 
                contact_type=contact_type),
            contact_constraint(get_body(mechanism, :sphere2), normal, 
                friction_coefficient=friction_coefficient, 
                contact_point=contact, 
                offset=[0.0, 0.0, radius * α], 
                contact_type=contact_type)
            ]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.0; 0.0; radius; zeros(3)])
        mechanism = Mechanism(origin, bodies, joints, contacts,
            gravity=gravity, 
            timestep=timestep)
    end

    return mechanism
end

function initialize_tippetop!(mechanism::Mechanism{T}; 
    x=zeros(3),
    q=one(UnitQuaternion), 
    v=zeros(3),
    ω=zeros(3)) where T

    joint2 = get_joint(mechanism, :fixed_joint)
    radius = joint2.constraints[1].vertices[1][3]
    origin = mechanism.origin
    pbody = get_body(mech, :sphere1)
    cbody = get_body(mech, :sphere2)

    zero_velocity!(mechanism)
   
    set_maximal_configurations!(origin, pbody; 
        parent_vertex=[0.0; 0.0; radius], 
        child_vertex=[0.0; 0.0; 0.0], 
        Δx=x, 
        Δq=q)
    set_maximal_configurations!(pbody,  cbody; 
        parent_vertex=[0.0; 0.0; radius], 
        child_vertex=[0.0; 0.0; 0.0], 
        Δx=[0.0; 0.0; 0.0], 
        Δq=one(UnitQuaternion))
    set_maximal_velocities!(origin, pbody; 
        parent_vertex=[0.0; 0.0; radius], 
        child_vertex=[0.0; 0.0; 0.0], 
        Δv=v, 
        Δω=ω)
    set_maximal_velocities!(pbody,  cbody; 
        parent_vertex=[0.0; 0.0; radius], 
        child_vertex=[0.0; 0.0; 0.0], 
        Δv=[0.0; 0.0; 0.0], 
        Δω=[0.0; 0.0; 0.0])
    return nothing
end
