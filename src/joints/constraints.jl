function constraint(mechanism, joint::JointConstraint)
    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)
    tra = constraint(joint.translational, pbody, cbody, joint.impulses[2][joint_impulse_index(joint,1)], mechanism.μ, mechanism.timestep)
    rot = constraint(joint.rotational, pbody, cbody, joint.impulses[2][joint_impulse_index(joint,2)], mechanism.μ, mechanism.timestep)
    return svcat(tra, rot)
end

function constraint_jacobian(joint::JointConstraint)
    tra = constraint_jacobian(joint.translational, joint.impulses[2][joint_impulse_index(joint, 1)])
    rot = constraint_jacobian(joint.rotational, joint.impulses[2][joint_impulse_index(joint, 2)])
    return diagonal_cat(tra, rot)
end

@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint, body::Body)
    relative = :(body.id == joint.parent_id ? :parent : :child)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(constraint_jacobian_configuration($relative,
        joint.translational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 1)], mechanism.timestep))
    rot = :(constraint_jacobian_configuration($relative,
        joint.rotational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 2)], mechanism.timestep))
    return :(vcat($tra, $rot))
end

# off-diagonal Jacobians
function off_diagonal_jacobians(mechanism, body::Body, joint::JointConstraint)
    return -impulse_map(mechanism, joint, body), constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep)
end

function off_diagonal_jacobians(mechanism, joint::JointConstraint, body::Body)
    return constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep), -impulse_map(mechanism, joint, body)
end

function off_diagonal_jacobians(mechanism, pbody::Body, cbody::Body)
    # time step
    timestep = mechanism.timestep

    # dimensions
    Ne = length(mechanism.joints)
    Nb = length(mechanism.bodies)
    Nc = length(mechanism.contacts)

    # Jacobian
    jacobian_parent_child = szeros(6, 6)
    jacobian_child_parent = szeros(6, 6)

    for connectionid in connections(mechanism.system, pbody.id)
        # joints
        if connectionid <= Ne
            joint = get_node(mechanism, connectionid)
            if pbody.id == joint.parent_id
                for element in (joint.translational, joint.rotational)
                    if cbody.id == joint.child_id
                        joint.damper && (jacobian_parent_child -= damper_jacobian_velocity(:parent, :child, element, pbody, cbody, timestep))
                        joint.damper && (jacobian_child_parent -= damper_jacobian_velocity(:child, :parent, element, pbody, cbody, timestep))
                    end
                end
            elseif cbody.id == joint.parent_id
                for element in (joint.translational, joint.rotational)
                    if pbody.id == joint.child_id
                        joint.damper && (jacobian_parent_child -= damper_jacobian_velocity(:parent, :child, element, cbody, pbody, timestep))
                        joint.damper && (jacobian_child_parent -= damper_jacobian_velocity(:child, :parent, element, cbody, pbody, timestep))
                    end
                end
            end
        end

        # contacts
        if connectionid > Ne + Nb
            contact = get_node(mechanism, connectionid)
            if pbody.id == contact.parent_id
                if cbody.id == contact.child_id
                    Jpc = impulse_map_jacobian(:parent, :child, contact.model,
                            pbody,
                            cbody,
                            contact.impulses[2],
                            mechanism.timestep) * integrator_jacobian_velocity(cbody, timestep)

                    Jcp = impulse_map_jacobian(:child, :parent, contact.model,
                            pbody,
                            cbody,
                            contact.impulses[2],
                            mechanism.timestep) * integrator_jacobian_velocity(pbody, timestep)
                    # impulse_map_jacobian_configuration(mechanism, body, contact) * integrator_jacobian_velocity(body, timestep)
                    # impulse_map(mechanism, contact, body) * contact.impulses[2]
                    jacobian_parent_child -= Jpc#damper_jacobian_velocity(:parent, :child, element, pbody, cbody, timestep)
                    jacobian_child_parent -= Jcp#damper_jacobian_velocity(:child, :parent, element, pbody, cbody, timestep)
                end
            elseif cbody.id == contact.parent_id
                if pbody.id == contact.child_id
                    Jpc = impulse_map_jacobian(:parent, :child, contact.model,
                            cbody,
                            pbody,
                            contact.impulses[2],
                            mechanism.timestep) * integrator_jacobian_velocity(pbody, timestep)

                    Jcp = impulse_map_jacobian(:child, :parent, contact.model,
                            cbody,
                            bbody,
                            contact.impulses[2],
                            mechanism.timestep) * integrator_jacobian_velocity(cbody, timestep)

                    jacobian_parent_child -= Jpc #damper_jacobian_velocity(:parent, :child, element, cbody, pbody, timestep)
                    jacobian_child_parent -= Jcp #damper_jacobian_velocity(:child, :parent, element, cbody, pbody, timestep)
                end
            end
        end
    end

    return jacobian_parent_child, jacobian_child_parent
end

# linear system
function set_matrix_vector_entries!(mechanism, matrix_entry::Entry, vector_entry::Entry, joint::JointConstraint)
    matrix_entry.value = constraint_jacobian(joint)
    vector_entry.value = -constraint(mechanism, joint)
end

function reset!(joint::JointConstraint{T,N,Nc}; scale::T=1.0) where {T,N,Nc}
    Nλ_tra = joint_length(joint.translational)
    Nb_tra = limits_length(joint.translational)
    Nλ_rot = joint_length(joint.rotational)
    Nb_rot = limits_length(joint.rotational)
    joint.impulses[1] = [scale * sones(2Nb_tra); szeros(Nλ_tra); scale * sones(2Nb_rot); szeros(Nλ_rot)]
    joint.impulses[2] = [scale * sones(2Nb_tra); szeros(Nλ_tra); scale * sones(2Nb_rot); szeros(Nλ_rot)]
    return
end
