function get_halfcheetah(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:halfcheetah,
    springs=0,
    dampers=0,
    parse_springs=true,
    parse_dampers=true,
    limits=true,
    joint_limits=Dict([
        (:bthigh, [-0.52,1.05]), 
        (:bshin, [-0.785,0.785]), 
        (:bfoot, [-0.400,0.785]), 
        (:fthigh, [-1,0.7]), 
        (:fshin, [-1.20,0.87]), 
        (:ffoot, [-0.5,0.5])]),
    keep_fixed_joints=false, 
    friction_coefficient=0.4,
    contact_feet=true,
    contact_body=true,
    T=Float64)
    
    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(String(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    # joint limits
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = Origin{T}()
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        normal = Z_AXIS
        names = contact_body ? getfield.(mechanism.bodies, :name) : [:ffoot, :bfoot]
        for name in names
            body = get_body(mechanism, name)
            if name == :torso # need special case for torso
                # torso
                pf = [+0.5 * body.shape.shapes[1].shapes[1].rh[2]; 0; 0]
                pb = [-0.5 * body.shape.shapes[1].shapes[1].rh[2]; 0; 0]
                o = body.shape.shapes[1].shapes[1].rh[1]
                push!(contacts, contact_constraint(body, normal; 
                    friction_coefficient, 
                    contact_origin=pf, 
                    contact_radius=o))
                push!(contacts, contact_constraint(body, normal; 
                    friction_coefficient, 
                    contact_origin=pb, 
                    contact_radius=o))

                # head
                pf = [+0.5 * body.shape.shapes[1].shapes[1].rh[2] + 0.214; 0; 0.1935]
                o = body.shape.shapes[1].shapes[1].rh[1]
                push!(contacts, contact_constraint(body, normal; 
                    friction_coefficient, 
                    contact_origin=pf, 
                    contact_radius=o))
            else
                p = [0;0; -0.5 * body.shape.shapes[1].rh[2]]
                o = body.shape.shapes[1].rh[1]
                push!(contacts, contact_constraint(body, normal;
                    friction_coefficient, 
                    contact_origin=p, 
                    contact_radius=o))
            end
        end
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0.576509, 0, 0.02792])

    # construction finished
    return mechanism
end

function initialize_halfcheetah!(mechanism::Mechanism{T}; 
    body_position=[0, 0],  
    body_orientation=0) where T

    set_minimal_coordinates!(mechanism,
                 get_joint(mechanism, :floating_joint),
                 [body_position[2] + 0.576509, -body_position[1], -body_orientation + 0.02792])
    for joint in mechanism.joints
        (joint.name != :floating_joint) && set_minimal_coordinates!(mechanism, joint, zeros(input_dimension(joint)))
    end
    zero_velocity!(mechanism)
end

function halfcheetahState(; x::T=0, z::T=0, θ::T=0) where T
    mechanism = get_mechanism(:halfcheetah)
    initialize!(mechanism, :halfcheetah, x=x, z=z, θ=θ)

    Nb = length(mechanism.bodies)
    x = zeros(13 * Nb)
    
    for (i, body) in enumerate(mechanism.bodies)
        x2 = body.state.x2
        v15 = zeros(3)
        q2 = body.state.q2
        ω15 = zeros(3)
        x[13 * (i-1) .+ (1:13)] = [x2;  v15; vector(q2); ω15]
    end
    return x
end