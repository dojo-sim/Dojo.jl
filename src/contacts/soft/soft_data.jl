
function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
        contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    Nd = data_dim(contact)
	∂imp∂θ = soft_impulse_jacobian_contact_data(mechanism, contact, body)
	∂dyn∂imp = -impulse_map(mechanism, contact, body)
	∂dyn∂θ = ∂dyn∂imp * ∂imp∂θ
	return ∂dyn∂θ
end

function contact_constraint_jacobian_contact_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs<:SoftContact{T,N}}
	soft_impulse_jacobian_contact_data(mechanism, contact, body)
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs}
	# [m; j; v15; ϕ15; x2; vector(q2)]
	∇ = szeros(T,6,13)
	∇xq = soft_impulse_jacobian_configuration(mechanism, contact, body)
	return [∇ ∇xq]
end
