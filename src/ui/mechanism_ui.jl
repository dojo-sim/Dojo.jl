function minimal_coordinates_velocities(mechanism::Mechanism)
    d = Dict()
    for joint in mechanism.joints
        push!(d, joint.id => [minimal_coordinates(mechanism, joint); minimal_velocities(mechanism, joint)])
    end
    return d
end

function minimal_coordinates(mechanism::Mechanism)
    d = Dict()
    for joint in mechanism.joints
        push!(d, joint.id => minimal_coordinates(mechanism, joint))
    end
    return d
end

function minimal_velocities(mechanism::Mechanism)
    d = Dict()
    for joint in mechanism.joints
        push!(d, joint.id => minimal_velocities(mechanism, joint))
    end
    return d
end

function minimal_configuration_vector(mechanism::Mechanism{T}) where T
    N = control_dimension(mechanism)
    x = zeros(T,N)
    off = 0
    for joint in mechanism.joints
        n = control_dimension(joint)
        x[off .+ (1:n)] += minimal_coordinates(mechanism, joint)
        off += n
    end
    return x
end

function minimal_velocity_vector(mechanism::Mechanism{T}) where T
    N = control_dimension(mech)
    x = zeros(T,N)
    off = 0
    for joint in mechanism.joints
        n = control_dimension(joint)
        x[off .+ (1:n)] += minimal_velocities(mechanism, joint)
        off += n
    end
    return x
end

function set_minimal_coordinates!(mechanism::Mechanism, dict)
    for (id,joint) in pairs(mechanism.joints)
        set_minimal_coordinates!(mechanism, joint, dict[id])
    end

    return
end

function set_minimal_velocities!(mechanism::Mechanism, dict)
    for (id,joint) in pairs(mechanism.joints)
        set_minimal_velocities!(mechanism, joint, dict[id])
    end

    return
end

function zero_velocity!(mechanism::Mechanism)
    # velocities
    for (i, body) in enumerate(mechanism.bodies)
        try
            set_maximal_velocities!(body, v=zeros(3), ω=zeros(3))
            set_previous_configuration!(body, mechanism.timestep)
        catch
            nothing
        end
    end
end

function set_input!(mechanism::Mechanism, dict)
    for (id,joint) in pairs(mechanism.joints)
        set_input!(joint, dict[id])
    end

    return
end

function set_current!(mechanism::Mechanism)
    for body in mechanism.bodies set_current!(body) end
end
