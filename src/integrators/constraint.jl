function constraint(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
    state = body.state
    timestep= mechanism.timestep

    mass = body.mass
    inertia = body.inertia
    gravity=mechanism.gravity

    x1, q1 = previous_configuration(state)
    x2, q2 = current_configuration(state)
    x3, q3 = next_configuration(state, timestep)

    # dynamics
    D1x = - 1.0 / timestep * mass * (x2 - x1) - 0.5 * timestep * mass * gravity
    D2x =   1.0 / timestep * mass * (x3 - x2) - 0.5 * timestep * mass * gravity
    D1q = -2.0 / timestep * LVᵀmat(q2)' * Lmat(q1) * Vᵀmat() * inertia * Vmat() * Lmat(q1)' * vector(q2)
    D2q = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)

    dynT = D2x + D1x
    dynR = D2q + D1q

    state.d = [dynT; dynR]

    # inputs
    state.d -= [state.F2; state.τ2]

    # impulses
    for id in connections(mechanism.system, body.id)
        Ne < id <= Ne + Nb && continue # body
        impulses!(mechanism, body, get_node(mechanism, id))
    end

    return state.d
end

# function constraint_jacobian_configuration(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}) where {T,Nn,Ne,Nb}
#     state = body.state
#     timestep = mechanism.timestep
#     mass = body.mass
#     inertia = body.inertia
#
#     # x1, q1 = previous_configuration(state)
#     x2, q2 = current_configuration(state)
#     x3, q3 = next_configuration(state, timestep)
#
#     # dynamics
#     dynT = I(3) * mass / timestep
#     dynR = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * (∂Rᵀmat∂q(Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)) + Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)')
#
#     Z33 = szeros(T, 3, 3)
#     Z34 = szeros(T, 3, 4)
#
#     state.D = [[dynT; Z33] [Z34; dynR]] * integrator_jacobian_velocity(body, timestep)
#     state.D += [[REG * I(3); Z33] [Z33; REG * I(3)]]
#
#     # inputs
#     nothing
#
#     # impulses
#     for id in connections(mechanism.system, body.id)
#         Ne < id <= Ne + Nb && continue # body
#         impulses_jacobian_velocity!(mechanism, body, get_node(mechanism, id))
#     end
#
#     return state.D
# end

function constraint_jacobian_configuration(mechanism::Mechanism{T,Nn,Ne,Nb}, body::Body{T}; reg::T=Dojo.REG) where {T,Nn,Ne,Nb}
    state = body.state
    timestep = mechanism.timestep
    mass = body.mass
    inertia = body.inertia

    x2, q2 = current_configuration(state)
    x3, q3 = next_configuration(state, timestep)

    I3 = SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    Z33 = szeros(T, 3, 3)
    Z34 = szeros(T, 3, 4)

    # dynamics
    dynT = I3 * mass / timestep
    dynR = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * (∂Rᵀmat∂q(Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)) + Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)')

    state.D = [[dynT; Z33] [Z34; dynR]] * integrator_jacobian_velocity(body, timestep)
    state.D += [[reg * I3; Z33] [Z33; reg * I3]]

    # inputs
    nothing

    # impulses
    for id in connections(mechanism.system, body.id)
        Ne < id <= Ne + Nb && continue # body
        impulses_jacobian_velocity!(mechanism, body, get_node(mechanism, id))
    end

    return state.D
end

function integrator_jacobian_velocity(body::Body{T}, timestep) where T
    state = body.state
    x2, v25, q2, ϕ25 = current_configuration_velocity(state)
    integrator_jacobian_velocity(x2, v25, q2, ϕ25, timestep)
end

function integrator_jacobian_configuration(body::Body{T},
        timestep; attjac::Bool=true) where T
    state = body.state
    x2, v25, q2, ϕ25 = current_configuration_velocity(state)
    integrator_jacobian_configuration(x2, v25, q2, ϕ25, timestep; attjac=attjac)
end

# linear system
function set_matrix_vector_entries!(mechanism, matrix_entry::Entry, vector_entry::Entry, body::Body)
    matrix_entry.value = constraint_jacobian_configuration(mechanism, body)
    vector_entry.value = -constraint(mechanism, body)
end
