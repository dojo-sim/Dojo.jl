function get_hopper(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:hopper,
    springs=10,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    limits=false,
    joint_limits=Dict([
        (:thigh, [0,150] * π / 180), 
        (:leg, [0,150] * π / 180), 
        (:foot, [-45,45] * π / 180)]),
    keep_fixed_joints=false, 
    friction_coefficient=2,
    contact_feet=true,
    contact_body=true,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]
    
    if contact_feet
        normal = Z_AXIS
        names = contact_body ? getfield.(mechanism.bodies, :name) : [:ffoot, :foot]
        for name in names
            body = get_body(mechanism, name)
            if name == :foot # need special case for foot
                # torso
                pf = [0, 0, +0.5 * body.shape.shapes[1].rh[2]]
                pb = [0, 0, -0.5 * body.shape.shapes[1].rh[2]]
                o = body.shape.shapes[1].rh[1]
                push!(contacts, contact_constraint(body, normal; 
                    friction_coefficient, 
                    contact_origin=pf, 
                    contact_radius=o))
                push!(contacts, contact_constraint(body, normal; 
                    friction_coefficient, 
                    contact_origin=pb, 
                    contact_radius=o))
            else
                p = [0; 0; 0.5 * body.shape.shapes[1].rh[2]]
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
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [1.25, 0, 0])

    # construction finished
    return mechanism
end

function initialize_hopper!(mechanism::Mechanism{T}; 
    body_position=[0, 0],
    body_orientation=0) where T
    #TODO add leg length

    set_minimal_coordinates!(mechanism,
                 get_joint(mechanism, :floating_joint),
                 [body_position[2] + 1.25 , -body_position[1], -body_orientation])
    for joint in mechanism.joints
        (joint.name != :floating_joint) && set_minimal_coordinates!(mechanism, joint, zeros(input_dimension(joint)))
    end
    zero_velocity!(mechanism)
end
