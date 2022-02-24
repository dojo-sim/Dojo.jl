get_joint_constraint(mechanism::Mechanism, id::Integer) = mechanism.joints[id]
get_body(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = id == 0 ? mechanism.origin : mechanism.bodies[id-Ne]
get_contact_constraint(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = mechanism.contacts[id-Ne-Nb]

function get_joint_constraint(mechanism::Mechanism, name::Symbol)
    for joint in mechanism.joints
        if joint.name == name
            return joint
        end
    end
    return
end

function get_body(mechanism::Mechanism, name::Symbol)
    if name == :origin
        return mechanism.origin
    else
        for body in mechanism.bodies
            if body.name == name
                return body
            end
        end
    end
    return
end

function get_contact_constraint(mechanism::Mechanism, name::Symbol)
    for contact in mechanism.contacts
        if contact.name == name
            return contact
        end
    end
    return
end

function get_node(mechanism::Mechanism{T,Nn,Ne,Nb}, id::Integer; origin::Bool=false) where {T,Nn,Ne,Nb}
    (origin && id == 0) && return mechanism.origin
    if id <= Ne
        return get_joint_constraint(mechanism, id)
    elseif id <= Ne+Nb
        return get_body(mechanism, id)
    else
        return get_contact_constraint(mechanism, id)
    end
end

function get_node(mechanism::Mechanism, name::Symbol)
    node = get_body(mechanism, name)
    if node === nothing
        node = get_joint_constraint(mechanism,name)
    end
    if node === nothing
        node = get_contact_constraint(mechanism,name)
    end
    return node
end

function initialize_state!(mechanism::Mechanism)
    for body in mechanism.bodies initialize_state!(body, mechanism.timestep) end
end

function off_diagonal_jacobians(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body1::Body, body2::Body) where {T,Nn,Ne,Nb,Ni}
    timestep = mechanism.timestep

    dimpulse_map_parentb = szeros(6, 6)
    dimpulse_map_childa = szeros(6, 6)

    for connectionid in connections(mechanism.system, body1.id)
        !(connectionid <= Ne) && continue # body
        joint = get_node(mechanism, connectionid)
        off = 0
        if body1.id == joint.parent_id
            for element in (joint.translational, joint.rotational)
                Nj = length(element)
                if body2.id == joint.child_id
                    joint.spring && (dimpulse_map_parentb -= spring_jacobian_velocity(:parent, :child, element, body1, body2, timestep)) #should be useless
                    joint.damper && (dimpulse_map_parentb -= damper_jacobian_velocity(:parent, :child, element, body1, body2, timestep))
                    joint.spring && (dimpulse_map_childa -= spring_jacobian_velocity(:child, :parent, element, body1, body2, timestep)) #should be useless
                    joint.damper && (dimpulse_map_childa -= damper_jacobian_velocity(:child, :parent, element, body1, body2, timestep))
                end
                off += Nj
            end
        elseif body2.id == joint.parent_id
            for element in (joint.translational, joint.rotational)
                Nj = length(element)
                if body1.id == joint.child_id
                    # joint.spring && (dimpulse_map_parentb -= spring_parent_jacobian_velocity_child(element, body2, body1, timestep)) #should be useless
                    joint.damper && (dimpulse_map_parentb -= damper_jacobian_velocity(:parent, :child, element, body2, body1, timestep))
                    # joint.spring && (dimpulse_map_childa -= spring_child_jacobian_velocity_parent(element, body2, body1, timestep)) #should be useless
                    joint.damper && (dimpulse_map_childa -= damper_jacobian_velocity(:child, :parent, element, body2, body1, timestep))
                end
                off += Nj
            end
        end
    end
    return dimpulse_map_parentb, dimpulse_map_childa
end
