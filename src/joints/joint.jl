################################################################################
# Joint
################################################################################
abstract type Joint{T,Nλ,Nb,N,Nb½} end

# joints
function joint_constraint(joint::Joint{T}, 
    xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T}
    return constraint_mask(joint) * displacement(joint, xa, qa, xb, qb)
end

function joint_constraint_jacobian_configuration(relative::Symbol, joint::Joint{T}, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, η) where {T}
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=false)
    return constraint_mask(joint) * [X Q]
end

# constraints
@inline function constraint(joint::Joint{T,Nλ,0}, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    joint_constraint(joint, xa, qa, xb, qb, η)
end

@inline constraint(joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)

# constraint Jacobians
@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ,0}, η) where {T,Nλ}
    return Diagonal(+1.00e-10 * sones(T,Nλ))
end

@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s, γ = split_impulses(joint, η)
    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline constraint_jacobian_configuration(relative::Symbol, joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint_jacobian_configuration(relative, joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)

@inline function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,0}, 
        xa::AbstractVector, qa::UnitQuaternion, 
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
end

@inline function constraint_jacobian_configuration(relative::Symbol, joint::Joint, body1::Node, body2::Node, child_id, λ, timestep)
    if body2.id == child_id
        return constraint_jacobian_configuration(relative, joint, 
                    next_configuration(body1.state, timestep)..., 
                    next_configuration(body2.state, timestep)..., 
                    λ)
    else
        return zero(joint)
    end
end

# masks
@inline constraint_mask(::Joint{T,0}) where T = szeros(T,0,3)
@inline constraint_mask(joint::Joint{T,1}) where T = joint.V3
@inline constraint_mask(joint::Joint{T,2}) where T = joint.V12
@inline constraint_mask(::Joint{T,3}) where T = SMatrix{3,3,T,9}(I)

@inline nullspace_mask(::Joint{T,0}) where T = SMatrix{3,3,T,9}(I)
@inline nullspace_mask(joint::Joint{T,1}) where T = joint.V12
@inline nullspace_mask(joint::Joint{T,2}) where T = joint.V3
@inline nullspace_mask(::Joint{T,3}) where T = szeros(T,0,3)

# impulse maps
@inline function impulse_map(relative::Symbol, joint::Joint, body1::Node, body2::Node, child_id, λ, timestep)
    if body2.id == child_id
        return impulse_map(relative, joint, current_configuration(body1.state)..., current_configuration(body2.state)..., λ)
    else
        return zero(joint)
    end
end

@inline function impulse_map(relative::Symbol, joint::Joint{T,Nλ,Nb}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb}
    return impulse_transform(relative, joint, xa, qa, xb, qb) * impulse_projector(joint)
end

@inline function impulse_map(relative::Symbol, joint::Joint{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    J = constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
    G = cat(Diagonal(sones(3)), LVᵀmat(relative == :parent ? qa : qb), dims=(1,2))
    return Diagonal([sones(3); 0.5 * sones(3)]) * transpose(J * G)
end

@inline function impulse_projector(joint::Joint{T,Nλ,0}) where {T,Nλ}
    zerodimstaticadjoint(constraint_mask(joint))
end

@inline function impulse_projector(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb}
    zerodimstaticadjoint([szeros(Nb,3); -nullspace_mask(joint); nullspace_mask(joint); constraint_mask(joint)])
end

# springs
@inline function spring_force(relative::Symbol, joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return spring_force(relative, joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

# dampers
@inline function damper_force(relative::Symbol, joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return damper_force(relative, joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

# inputs
@inline function set_input!(joint::Joint, input::SVector)
    joint.input = zerodimstaticadjoint(nullspace_mask(joint)) * input
    return
end

@inline set_input!(joint::Joint) = return

@inline function add_input!(joint::Joint, input::SVector)
    joint.input += zerodimstaticadjoint(nullspace_mask(joint)) * input
    return
end

@inline add_input!(joint::Joint) = return

@inline function input_jacobian_control(relative::Symbol, joint::Joint, body1::Node, body2::Node, child_id)
    return input_jacobian_control(relative, joint, current_configuration(body1.state)..., current_configuration(body2.state)...) * zerodimstaticadjoint(nullspace_mask(joint))
end

# minimal coordinates
@inline minimal_coordinates(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 3 - Nλ)

################################################################################
# Utilities
################################################################################
Base.length(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
Base.zero(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, Nλ, 6)
joint_length(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
limits_length(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb} = Nb
impulses_length(joint::Joint{T,Nλ,Nb,N}) where {T,Nλ,Nb,N} = N

function split_impulses(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s = η[SVector{Nb,Int}(1:Nb)]
    γ = η[SVector{Nb,Int}(Nb .+ (1:Nb))]
    return s, γ
end

function joint_impulse_index(joint::Joint{T,Nλ,Nb,N}, s::Int) where {T,Nλ,Nb,N}
    ind = SVector{N,Int}(s+1:s+N)
    return ind
end