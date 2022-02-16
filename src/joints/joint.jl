abstract type Joint{T,Nλ,Nb,N,Nb½} end

getT(joint::Joint{T}) where T = T
Base.length(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
Base.zero(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, Nλ, 6)

λlength(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
blength(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb} = Nb
ηlength(joint::Joint{T,Nλ,Nb,N}) where {T,Nλ,Nb,N} = N

function get_sγ(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s = η[SVector{Nb,Int}(1:Nb)]
    γ = η[SVector{Nb,Int}(Nb .+ (1:Nb))]
    return s, γ
end

function λindex(joint::Joint{T,Nλ,Nb,N}, s::Int) where {T,Nλ,Nb,N}
    ind = SVector{N,Int}(s+1:s+N)
    return ind
end

# Joint constraints
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

@inline function constraint(joint::Joint{T,Nλ,0}, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    joint_constraint(joint, xa, qa, xb, qb, η)
end

@inline function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,0}, 
        xa::AbstractVector, qa::UnitQuaternion, 
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
end

@inline function constraint_jacobian_configuration(relative::Symbol, joint::Joint, body1::Node, body2::Node, child_id, λ, timestep)
    if body2.id == child_id
        return constraint_jacobian_configuration(relative, joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)
    else
        return zero(joint)
    end
end

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

# With joint limits
@inline function impulse_projector(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb}
    zerodimstaticadjoint([szeros(Nb,3); -nullspace_mask(joint); nullspace_mask(joint); constraint_mask(joint)])
end

# Without joint limits
@inline function impulse_projector(joint::Joint{T,Nλ,0}) where {T,Nλ}
    zerodimstaticadjoint(constraint_mask(joint))
end

@inline function spring_parent(joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return spring_parent(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function spring_child(joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return spring_child(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function damper_parent(joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return damper_parent(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function damper_child(joint::Joint, body1::Node, body2::Node, timestep, child_id; unitary::Bool=false)
    if body2.id == child_id
        return damper_child(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function apply_input!(joint::Joint, body1::Node, body2::Node, timestep::T, clear::Bool) where T
    apply_input!(joint, body1.state, body2.state, timestep, clear)
    return
end

@inline constraint_mask(::Joint{T,0}) where T = szeros(T,0,3)
@inline nullspace_mask(::Joint{T,0}) where T = SMatrix{3,3,T,9}(I)
@inline constraint_mask(joint::Joint{T,1}) where T = joint.V3
@inline nullspace_mask(joint::Joint{T,1}) where T = joint.V12
@inline constraint_mask(joint::Joint{T,2}) where T = joint.V12
@inline nullspace_mask(joint::Joint{T,2}) where T = joint.V3
@inline constraint_mask(::Joint{T,3}) where T = SMatrix{3,3,T,9}(I)
@inline nullspace_mask(::Joint{T,3}) where T = szeros(T,0,3)

@inline constraint(joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)

@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ,0}, η) where {T,Nλ}
    return Diagonal(+1.00e-10 * sones(T,Nλ))
end

@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s, γ = get_sγ(joint, η)
    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline constraint_jacobian_configuration(relative::Symbol, joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint_jacobian_configuration(relative, joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)


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

@inline function input_jacobian_control_parent(joint::Joint, body1::Node, body2::Node, child_id)
    return input_jacobian_control_parent(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspace_mask(joint))
end

@inline function input_jacobian_control_child(joint::Joint{T,Nλ}, body1::Node, body2::Node, child_id) where {T,Nλ}
    if body2.id == child_id
        return input_jacobian_control_child(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspace_mask(joint))
    else
        return szeros(T, 6, 3 - Nλ)
    end
end

@inline minimal_coordinates(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 3 - Nλ)