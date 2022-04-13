# previous
previous_configuration(state::State) = (state.x1, state.q1)
previous_configuration_velocity(state::State) = (state.x1, state.v15, state.q1, state.ϕ15)

# current
current_position(state::State) = state.x2
current_orientation(state::State) = state.q2
current_configuration(state::State) = (current_position(state), current_orientation(state))
current_velocity(state::State) = (state.vsol[2], state.ϕsol[2])
current_configuration_velocity(state::State) = (current_position(state), state.vsol[2], current_orientation(state), state.ϕsol[2])
initial_configuration_velocity(state::State) = (current_position(state), state.v15, current_orientation(state), state.ϕ15)

# next
next_position(x2::SVector{3}, v25::SVector{3}, timestep::Real) = x2 + v25 * timestep
next_orientation(q2::Quaternion, ϕ25::SVector{3}, timestep::Real) = q2 * quaternion_map(ϕ25, timestep) * timestep / 2
next_position(state::State, timestep) = next_position(state.x2, state.vsol[2], timestep)
next_orientation(state::State, timestep) = next_orientation(state.q2, state.ϕsol[2], timestep)
next_configuration(state::State, timestep) = (next_position(state, timestep), next_orientation(state, timestep))
next_configuration_velocity(state::State, timestep) = (next_position(state, timestep), state.vsol[2], next_orientation(state, timestep), state.ϕsol[2])

# angular velocity
function angular_velocity(q1::Quaternion, q2::Quaternion, timestep)
    2.0 / timestep  * Vmat() * Lᵀmat(q1) * vector(q2)
end

function ∂angular_velocity∂q1(q1::Quaternion, q2::Quaternion, timestep)
    2.0 / timestep  * Vmat() * Rmat(q2) * Tmat()
end

function ∂angular_velocity∂q2(q1::Quaternion, q2::Quaternion, timestep)
    2.0 / timestep  * Vmat() * Lᵀmat(q1)
end

# Jacobians
function integrator_jacobian_velocity(x2::AbstractVector, v25::AbstractVector, q2::Quaternion, ϕ25::SVector{3}, timestep::T) where T
    V = [linear_integrator_jacobian_velocity(x2, v25, timestep) szeros(T,3,3)]
    Ω = [szeros(T,4,3) rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)]
    return [V; Ω] # 7x6
end

function integrator_jacobian_configuration(x2::AbstractVector, v25::AbstractVector, q2::Quaternion, ϕ25::SVector{3}, timestep::T;
    attjac::Bool=true) where T

    Z = attjac ? szeros(T,3,3) : szeros(T,3,4)
    X = [linear_integrator_jacobian_position(x2, v25, timestep) Z]
    Q = [szeros(T,4,3) rotational_integrator_jacobian_orientation(q2, ϕ25, timestep; attjac=attjac)]
    return [X; Q] # 7x6 or 7x7
end

function linear_integrator_jacobian_position(x2::AbstractVector, v2::AbstractVector, timestep::T) where T
    return SMatrix{3,3,T,9}(Diagonal(sones(T, 3)))
end

function linear_integrator_jacobian_velocity(x2::AbstractVector, v2::AbstractVector, timestep::T) where T
    return timestep * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
end

function rotational_integrator_jacobian_orientation(q2::Quaternion, ϕ25::SVector{3}, timestep::Real;
    attjac::Bool = true)
    M = Rmat(quaternion_map(ϕ25, timestep) * timestep / 2)
    attjac && (M *= LVᵀmat(q2))
    return M
end

function rotational_integrator_jacobian_velocity(q2::Quaternion, ϕ25::SVector{3}, timestep::Real)
    return Lmat(q2) * quaternion_map_jacobian(ϕ25, timestep) * timestep / 2
end
