function correction!(mechanism)
	system = mechanism.system
	residual_entries = mechanism.residual_entries
    for id in reverse(system.dfs_list)
        node = get_node(mechanism, id)
        correction!(mechanism, residual_entries[id], get_entry(system, id), node)
    end
	return
end

correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, node::Node) = nothing

function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, ::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	μ = mechanism.μ
	residual_entry.value += [- Δs .* Δγ .+ μ; szeros(N½)]
    return
end

function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, contact::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
	cont = contact.model
	μ = mechanism.μ
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	residual_entry.value += [[-Δs[1] * Δγ[1]; -cone_product(Δs[2:4], Δγ[2:4])] + μ * neutral_vector(cont); szeros(N½)]
    return
end

function correction!(mechanism::Mechanism{T}, residual_entry::Entry, step_entry::Entry, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
	cor = correction(mechanism, step_entry, joint)
	residual_entry.value += cor
    return
end

@generated function correction(mechanism::Mechanism{T}, step_entry::Entry, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    cor_tra = :(correction(joint.translational, step_entry.value[joint_impulse_index(joint, 1)], mechanism.μ))
    cor_rot = :(correction(joint.rotational, step_entry.value[joint_impulse_index(joint, 2)], mechanism.μ))
    return :(vcat($cor_tra, $cor_rot))
end

function correction(joint::Joint{T,Nλ,Nb,N}, Δ, μ) where {T,Nλ,Nb,N}
    Δs, Δγ = split_impulses(joint, Δ)
	return [- Δs .* Δγ .+ μ; szeros(Nb + Nλ)]
end