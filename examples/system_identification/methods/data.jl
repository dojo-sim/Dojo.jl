# simdata_dim(mechanism::Mechanism) = sum(simdata_dim.(mech.contacts))
# simdata_dim(contact::ContactConstraint) = sum(simdata_dim.(contact.constraints))
# simdata_dim(model::NonlinearContact) = 5 # [friction_coefficient, p, contact_radius]
# simdata_dim(model::LinearContact) = 5 # [friction_coefficient, p, contact_radius]
# simdata_dim(model::ImpactContact) = 4 # [p, contact_radius]
#
# function set_simulator_data!(model::NonlinearContact, data::AbstractVector)
# 	model.friction_coefficient = data[1]
#     model.collision.contact_radius = norm(data[SVector{3,Int}(2:4)]) #TODO: Confirm this is correct
#     model.collision.contact_origin = data[SVector{3,Int}(5:7)]
#     return nothing
# end
#
# function set_simulator_data!(model::LinearContact, data::AbstractVector)
# 	model.friction_coefficient = data[1]
#     model.collision.contact_radius = norm(data[SVector{3,Int}(2:4)]) #TODO: Confirm this is correct
#     model.collision.contact_origin = data[SVector{3,Int}(5:7)]
#     return nothing
# end
#
# function set_simulator_data!(model::ImpactContact, data::AbstractVector)
#     model.collision.contact_radius = norm(data[SVector{3,Int}(2:4)]) #TODO: Confirm this is correct
#     model.collision.contact_origin = data[SVector{3,Int}(5:7)]
#     return nothing
# end
#
# function set_simulator_data!(mechanism::Mechanism, data::AbstractVector)
#     c = 0
#     for contact in mechanism.contacts
# 		for element in contact.constraints
# 			N = simdata_dim(element)
# 	        set_simulator_data!(element, data[c .+ (1:N)]); c += N
# 		end
#     end
#     return nothing
# end
#
# function get_simulator_data(model::NonlinearContact)
# 	return [model.friction_coefficient; model.collision.contact_radius; model.collision.contact_origin]
# end
#
# function get_simulator_data(modelset::LinearContact)
# 	return [model.friction_coefficient; model.collision.contact_radius; model.collision.contact_origin]
# end
#
# function get_simulator_data(model::ImpactContact)
# 	return [model.collision.contact_radius; model.collision.contact_origin]
# end
#
# function get_simulator_data(mechanism::Mechanism{T}) where T
#     data = zeros(T,simdata_dim(mechanism))
# 	e = 0
#     for contact in mechanism.contacts
# 		for element in contact.constraints
# 			N = simdata_dim(element)
#         	data[e .+ (1:N)] = get_simulator_data(model); e += N
# 		end
#     end
#     return data
# end
