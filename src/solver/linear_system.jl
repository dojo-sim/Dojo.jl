function set_entries!(mechanism::Mechanism)
    system = mechanism.system

    for id in reverse(system.dfs_list)
        for child_id in system.cyclic_children[id]
            zero_LU!(get_entry(system, id, child_id), get_entry(system, child_id, id))
        end

        node = get_node(mechanism, id)
        set_matrix_vector_entries!(mechanism, get_entry(system, id, id), get_entry(system, id), node)

        for child_id in children(system,id)
            set_LU!(mechanism, get_entry(system, id, child_id), get_entry(system, child_id, id), node, get_node(mechanism, child_id))
        end
    end
    return
end

function set_LU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, nodea::Node, nodeb::Node)
    L, U = off_diagonal_jacobians(mechanism, nodea, nodeb)
    matrix_entry_L.value = L
    matrix_entry_U.value = U
    return
end

function zero_LU!(matrix_entry_L::Entry, matrix_entry_U::Entry)
    matrix_entry_L.value *= 0
    matrix_entry_U.value *= 0
    return
end

function pull_residual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.residual_entries[i].value = mechanism.system.vector_entries[i].value
	end
	return
end

function push_residual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.system.vector_entries[i].value = mechanism.residual_entries[i].value
	end
	return
end

function pull_matrix!(mechanism::Mechanism)
	mechanism.matrix_entries.nzval .= mechanism.system.matrix_entries.nzval #TODO: make allocation free
	return
end

function push_matrix!(mechanism::Mechanism)
	mechanism.system.matrix_entries.nzval .= mechanism.matrix_entries.nzval #TODO: make allocation free
	return
end

function update!(body::Body)
    body.state.vsol[1] = body.state.vsol[2]
    body.state.ωsol[1] = body.state.ωsol[2]
    return
end

function update!(joint::JointConstraint)
    joint.impulses[1] = joint.impulses[2]
    return
end

function update!(contact::ContactConstraint)
    contact.impulses_dual[1] = contact.impulses_dual[2]
    contact.impulses[1] = contact.impulses[2]
    return
end


