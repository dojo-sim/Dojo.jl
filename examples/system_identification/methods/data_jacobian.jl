#
# function simdata_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
# 	nsol = solution_dimension(mechanism)
#     ndata = simdata_dim(mechanism)
#     data = get_simulator_data(mechanism)
#
# 	A = zeros(nsol, ndata)
# 	njoints = joint_dimension(mechanism)
# 	n = sum(length.(mech.joints.values))
# 	rb = njoints # row contact
# 	ri = njoints + 6Nb # row body
# 	c = 0 # column
# 	for contact in mechanism.contacts
# 		for element in contact.constraints
# 			N = length(contact)
# 			N½ = Int(N/2)
# 			Ns = simdata_dim(element)
# 			∇contact, ∇body = ∂g∂simdata(mechanism, contact) # assumes only one element per contact
# 			A[rb + 3 .+ (1:3), c .+ (1:Ns)] -= ∇body # only works for one body for now
# 			A[ri + N½ .+ (1:N½), c .+ (1:Ns)] -= ∇contact
# 			ri += N
# 			c += Ns
# 		end
# 	end
# 	return A
# end
#
# function ∂g∂simdata(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:NonlinearContact{T,N}}
#     model = contact.model
# 	p = model.collision.contact_origin
# 	offset = model.collision.contact_normal' * model.collision.contact_radius
#     body = get_body(mechanism, contact.parent_id)
#     x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
#     x3, q3 = next_configuration(body.state, mechanism.timestep)
# 	s = contact.impulses_dual[2]
# 	γ = contact.impulses[2]
#
# 	# Contribution to Injointonstraint
# 	∇friction_coefficient = SA[0,γ[1],0,0]
# 	∇off = [-model.collision.contact_normal; szeros(T,1,3); -model.contact_tangent * skew(vector_rotate(ϕ25, q3))]
# 	∇p = [model.collision.contact_normal * rotation_matrix(q3); szeros(T,1,3); model.contact_tangent * skew(vector_rotate(ϕ25, q3)) * rotation_matrix(q3)]
# 	∇contact = [∇friction_coefficient ∇off ∇p]
#
# 	# Contribution to Body dynamics
# 	∇friction_coefficient = szeros(T,3)
# 	X = force_mapping(model, x3, q3)
# 	# this what we differentiate: Qᵀγ = - skew(p - vector_rotate(offset, inv(q3))) * VRmat(q3) * LᵀVᵀmat(q3) * X' * γ
# 	∇off = - ∂skew∂p(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ) * -rotation_matrix(inv(q3))
# 	∇p = - ∂skew∂p(VRmat(q3) * LᵀVᵀmat(q3) * X' * γ)
# 	∇body = [∇friction_coefficient ∇off ∇p]
# 	return ∇contact, ∇body
# end
