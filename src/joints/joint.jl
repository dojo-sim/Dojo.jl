abstract type Joint{T,Nλ,Nb,N} end

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

@inline function impulse_map_parent(joint::Joint, body1::Node, body2::Node, childid, λ, timestep)
    if body2.id == childid
        return impulse_map_parent(joint, current_configuration(body1.state)..., current_configuration(body2.state)..., λ)
    else
        return zero(joint)
    end
end

@inline function impulse_map_child(joint::Joint, body1::Node, body2::Node, childid, λ, timestep)
    if body2.id == childid
        return impulse_map_child(joint, current_configuration(body1.state)..., current_configuration(body2.state)..., λ)
    else
        return zero(joint)
    end
end

function impulse_map_parent(joint::Joint, statea::State, stateb::State, η, timestep)
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_parent(joint, xa, qa, xb, qb, η)
end

function impulse_map_child(joint::Joint, statea::State, stateb::State, η, timestep)
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_child(joint, xa, qa, xb, qb, η)
end

@inline function impulse_map_parent(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    return impulse_transform_parent(joint, xa, qa, xb, qb) * impulse_projector(joint)
end

@inline function impulse_map_child(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    return impulse_transform_child(joint, xa, qa, xb, qb) * impulse_projector(joint)
end

# With joint limits
@inline function impulse_projector(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb}
    zerodimstaticadjoint([szeros(Nb,3); -nullspace_mask(joint); nullspace_mask(joint); constraint_mask(joint)])
end

# Without joint limits
@inline function impulse_projector(joint::Joint{T,Nλ,0}) where {T,Nλ}
    zerodimstaticadjoint(constraint_mask(joint))
end


################################################################################
# Derivatives
################################################################################

function impulse_map_parent_jacobian_parent(joint::Joint, pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(Ga*λ)/∂(xa,qa)
    p = impulse_projector(joint) * λ
    impulse_transform_parent_jacobian_parent(joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end

function impulse_map_parent_jacobian_child(joint::Joint, pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(Ga*λ)/∂(xb,qb)
    p = impulse_projector(joint) * λ
    impulse_transform_parent_jacobian_child(joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end

function impulse_map_child_jacobian_parent(joint::Joint, pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(Gb*λ)/∂(xa,qa)
    p = impulse_projector(joint) * λ
    impulse_transform_child_jacobian_parent(joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end

function impulse_map_child_jacobian_child(joint::Joint, pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(Gb*λ)/∂(xb,qb)
    p = impulse_projector(joint) * λ
    impulse_transform_child_jacobian_child(joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end

## Discrete-time velocity derivatives (for dynamics)
@inline function constraint_jacobian_parent(joint::Joint, body1::Node, body2::Node, childid, λ, timestep)
    if body2.id == childid
        return constraint_jacobian_parent(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)
    else
        return zero(joint)
    end
end

@inline function constraint_jacobian_child(joint::Joint, body1::Node, body2::Node, childid, λ, timestep)
    if body2.id == childid
        return constraint_jacobian_child(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)

    else
        return zero(joint)
    end
end

@inline function spring_parent(joint::Joint, body1::Node, body2::Node, timestep, childid; unitary::Bool=false)
    if body2.id == childid
        return spring_parent(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function spring_child(joint::Joint, body1::Node, body2::Node, timestep, childid; unitary::Bool=false)
    if body2.id == childid
        return spring_child(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function damper_parent(joint::Joint, body1::Node, body2::Node, timestep, childid; unitary::Bool=false)
    if body2.id == childid
        return damper_parent(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function damper_child(joint::Joint, body1::Node, body2::Node, timestep, childid; unitary::Bool=false)
    if body2.id == childid
        return damper_child(joint, body1, body2, timestep, unitary=unitary)
    else
        return szeros(T, 6)
    end
end

@inline function apply_input!(joint::Joint, body1::Node, body2::Node, timestep::T, clear::Bool) where T
    apply_input!(joint, body1.state, body2.state, timestep, clear)
    return
end

Joint0 = Joint{T,0} where T
Joint1 = Joint{T,1} where T
Joint2 = Joint{T,2} where T
Joint3 = Joint{T,3} where T

@inline constraint_mask(::Joint0{T}) where T = szeros(T,0,3)
@inline nullspace_mask(::Joint0{T}) where T = SMatrix{3,3,T,9}(I)
@inline constraint_mask(joint::Joint1) = joint.V3
@inline nullspace_mask(joint::Joint1) = joint.V12
@inline constraint_mask(joint::Joint2) = joint.V12
@inline nullspace_mask(joint::Joint2) = joint.V3
@inline constraint_mask(::Joint3{T}) where T = SMatrix{3,3,T,9}(I)
@inline nullspace_mask(::Joint3{T}) where T = szeros(T,0,3)

@inline constraint(joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)

@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ}, λ) where {T,Nλ}
    return Diagonal(+1.00e-10 * sones(T,Nλ))
end

@inline constraint_jacobian_parent(joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint_jacobian_parent(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)
@inline constraint_jacobian_child(joint::Joint, body1::Node, body2::Node, λ, timestep) = constraint_jacobian_child(joint, next_configuration(body1.state, timestep)..., next_configuration(body2.state, timestep)..., λ)


@inline function set_input!(joint::Joint, Fτ::SVector)
    joint.Fτ = zerodimstaticadjoint(nullspace_mask(joint)) * Fτ
    return
end

@inline set_input!(joint::Joint) = return

@inline function add_force!(joint::Joint, Fτ::SVector)
    joint.Fτ += zerodimstaticadjoint(nullspace_mask(joint)) * Fτ
    return
end

@inline add_force!(joint::Joint) = return

@inline function input_jacobian_control_parent(joint::Joint, body1::Node, body2::Node, timestep, childid)
    return input_jacobian_control_parent(joint, body1.state, body2.state, timestep) * zerodimstaticadjoint(nullspace_mask(joint))
end

@inline function input_jacobian_control_child(joint::Joint{T,Nλ}, body1::Node, body2::Node, timestep, childid) where {T,Nλ}
    if body2.id == childid
        return input_jacobian_control_child(joint, body1.state, body2.state, timestep) * zerodimstaticadjoint(nullspace_mask(joint))
    else
        return szeros(T, 6, 3 - Nλ)
    end
end

@inline minimal_coordinates(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 3 - Nλ)

function add_limits(mech::Mechanism, eq::JointConstraint;
    # NOTE: this only works for joints between serial chains (ie, single child joints)
    tra_limits=eq.constraints[1].joint_limits,
    rot_limits=eq.constraints[1].joint_limits)

    # update translational
    tra = eq.constraints[1]
    T = typeof(tra).parameters[1]
    Nλ = typeof(tra).parameters[2]
    Nb½ = length(tra_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    tra_limit = (Translational{T,Nλ,Nb,N,Nb½,N̄λ}(tra.V3, tra.V12, tra.vertices, tra.spring, tra.damper, tra.spring_offset, tra_limits, tra.spring_type, tra.Fτ), eq.parentid, eq.childids[1])

    # update rotational
    rot = eq.constraints[2]
    T = typeof(rot).parameters[1]
    Nλ = typeof(rot).parameters[2]
    Nb½ = length(rot_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    rot_limit = (Rotational{T,Nλ,Nb,N,Nb½,N̄λ}(rot.V3, rot.V12, rot.qoffset, rot.spring, rot.damper, rot.spring_offset, rot_limits, rot.spring_type, rot.Fτ), eq.parentid, eq.childids[1])
    JointConstraint((tra_limit, rot_limit); name=eq.name)
end
