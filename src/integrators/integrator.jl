@inline previous_configuration(state::State) = (state.x1, state.q1)
@inline previous_configuration_velocity(state::State) = (state.x1, state.v15, state.q1, state.ϕ15)

# Initial conditions of a body
@inline initial_configuration_velocity(state::State) = (current_position(state, k=1), state.v15, current_orientation(state, k=1), state.ϕ15)

@inline current_position(state::State; k=1) = state.x2[k]
@inline current_orientation(state::State; k=1) = state.q2[k]
@inline current_configuration(state::State; k=1) = (current_position(state, k=k), current_orientation(state, k=k))
@inline current_velocity(state::State) = (state.vsol[2], state.ϕsol[2])
@inline current_configuration_velocity(state::State) = (current_position(state, k=1), state.vsol[2], current_orientation(state, k=1), state.ϕsol[2])

@inline next_position(x2::SVector{3,T}, v25::SVector{3,T}, timestep::T) where T = x2 + v25 * timestep
@inline next_orientation(q2::UnitQuaternion{T}, ϕ25::SVector{3,T}, timestep::T) where T = q2 * quaternion_map(ϕ25, timestep) * timestep / 2
@inline next_position(state::State, timestep) = next_position(state.x2[1], state.vsol[2], timestep)
@inline next_orientation(state::State, timestep) = next_orientation(state.q2[1], state.ϕsol[2], timestep)
@inline next_configuration(state::State, timestep) = (next_position(state, timestep), next_orientation(state, timestep))

@inline function quaternion_map_jacobian(ω::SVector{3}, timestep)
    msq = -sqrt(4 / timestep^2 - dot(ω, ω))
    return [ω' / msq; I]
end

@inline function quaternion_map(ω, timestep)
    return UnitQuaternion(sqrt(4 / timestep^2 - dot(ω, ω)), ω, false)
end

function cayley(ω)
    UnitQuaternion(1.0 / sqrt(1.0 + norm(ω)^2.0) * [1.0; ω], false)
end

function cayley_jacobian(ω)
    ω₁, ω₂, ω₃ = ω
    a = sqrt(1.0 + sqrt(abs2(ω₁) + abs2(ω₂) + abs2(ω₃))^2.0)^-3
    b = sqrt(1.0 + sqrt(abs2(ω₁) + abs2(ω₂) + abs2(ω₃))^2.0)^-1
    SMatrix{4,3}([
                 -ω₁*a -ω₂*a -ω₃*a;
                 (b - (ω₁^2)*a) (-ω₁ * ω₂ * a) (-ω₁ * ω₃ * a);
                 (-ω₁ * ω₂ * a) (b - (ω₂^2)*a) (-ω₂ * ω₃ * a);
                 (-ω₁ * ω₃ * a) (-ω₂ * ω₃ * a) (b - (ω₃^2)*a);
                 ])
end

# I think this is the inverse of next_orientation, we recover ϕ15 from q1, q2 and h
function angular_velocity(q1::UnitQuaternion, q2::UnitQuaternion, timestep)
    2.0 / timestep  * Vmat() * Lᵀmat(q1) * vector(q2)
end

@inline function initialize_state!(body::Body{T}, timestep) where T
    state = body.state
    x2 = state.x2[1]
    q2 = state.q2[1]
    v15 = state.v15
    ϕ15 = state.ϕ15

    state.x1 = x2 - v15*timestep
    state.q1 = q2 * quaternion_map(-ϕ15,timestep) * timestep / 2

    state.F2[1] = szeros(T,3)
    state.τ2[1] = szeros(T,3)

    return
end

@inline function update_state!(body::Body{T}, timestep) where T
    state = body.state

    state.x1 = state.x2[1]
    state.q1 = state.q2[1]

    state.v15 = state.vsol[2]
    state.ϕ15 = state.ϕsol[2]

    state.x2[1] = state.x2[1] + state.vsol[2]*timestep
    state.q2[1] = state.q2[1] * quaternion_map(state.ϕsol[2], timestep) * timestep / 2

    state.F2[1] = szeros(T,3)
    state.τ2[1] = szeros(T,3)
    return
end

@inline function set_solution!(body::Body)
    state = body.state
    state.vsol[1] = state.v15
    state.vsol[2] = state.v15
    state.ϕsol[1] = state.ϕ15
    state.ϕsol[2] = state.ϕ15
    return
end

function integrator_jacobian_velocity(q2::UnitQuaternion{T}, ϕ25::SVector{3,T}, timestep::T) where T
    V = [linear_integrator_jacobian_velocity(timestep) szeros(T,3,3)]
    Ω = [szeros(T,4,3) rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)]
    return [V; Ω] # 7x6
end

function integrator_jacobian_configuration(q2::UnitQuaternion{T}, ϕ25::SVector{3,T},
        timestep::T; attjac::Bool=true) where T
    Z = attjac ? szeros(T,3,3) : szeros(T,3,4)
    X = [linear_integrator_jacobian_position() Z]
    Q = [szeros(T,4,3) rotational_integrator_jacobian_orientation(q2, ϕ25, timestep; attjac=attjac)]
    return [X; Q] # 7x6 or 7x7
end

function linear_integrator_jacobian_position(; T=Float64)
    return SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
end

function linear_integrator_jacobian_velocity(timestep::T) where T
    return timestep * SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
end

function rotational_integrator_jacobian_orientation(q2::UnitQuaternion{T}, ϕ25::SVector{3,T}, timestep::T; attjac::Bool = true) where T
    M = Rmat(quaternion_map(ϕ25, timestep) * timestep/2)
    attjac && (M *= LVᵀmat(q2))
    return M
end

function rotational_integrator_jacobian_velocity(q2::UnitQuaternion{T}, ϕ25::SVector{3,T}, timestep::T) where T
    return Lmat(q2) * quaternion_map_jacobian(ϕ25, timestep) * timestep/2
end
