
function body_constraint_jacobian_contact_data(mechanism::Mechanism, body::Body{T},
        contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    Nd = data_dim(contact)
	∂imp∂θ = soft_impulse_jacobian_contact_data(mechanism, contact, body)
	∂dyn∂imp = -impulse_map(mechanism, contact, body)
	∂dyn∂θ = ∂dyn∂imp * ∂imp∂θ
	# @warn "setting this to zero"
	return ∂dyn∂θ * 1.00000000000
end

function contact_constraint_jacobian_contact_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs<:SoftContact{T,N}}
	# @warn "setting this to zero"
	soft_impulse_jacobian_contact_data(mechanism, contact, body) * 1.000000000000
end

function contact_constraint_jacobian_body_data(mechanism::Mechanism,
		contact::SoftContactConstraint{T,N,Nc,Cs}, body::Body{T}) where {T,N,Nc,Cs}
	# [m; j; v15; ϕ15; x2; vector(q2)]
	# soft_impulse only depends on x2, q2, v25, ϕ25 and not v15, ϕ15
	∇ = szeros(T,6,13)
	# ∇vϕ = soft_impulse_jacobian_velocity(mechanism, contact, body)
	∇xq = soft_impulse_jacobian_configuration(mechanism, contact, body)
	# return [∇ ∇vϕ ∇xq]
	return [∇ ∇xq]
end
