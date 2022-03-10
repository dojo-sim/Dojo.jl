"""
   Joint{T} 

   Abstract type for 3-dimensional constraint between two Body objects
"""
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
function constraint(joint::Joint{T,Nλ,0},
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion, η, μ) where {T,Nλ}
    joint_constraint(joint, xa, qa, xb, qb, η)
end

constraint(joint::Joint, pbody::Node, cbody::Node, λ, μ, timestep) = constraint(joint, next_configuration(pbody.state, timestep)..., next_configuration(cbody.state, timestep)..., λ, μ)

# constraint Jacobians
function constraint_jacobian(joint::Joint{T,Nλ,0}, η) where {T,Nλ}
    return Diagonal(REG * sones(T,Nλ))
end

function constraint_jacobian(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s, γ = split_impulses(joint, η)
    c1 = [Diagonal(γ + REG * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + REG * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(REG * sones(T, Nλ))]
    return [c1 c2 c3]
end

constraint_jacobian_configuration(relative::Symbol, joint::Joint, pbody::Node, cbody::Node, λ, timestep) = constraint_jacobian_configuration(relative, joint, next_configuration(pbody.state, timestep)..., next_configuration(cbody.state, timestep)..., λ)

function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
end

# masks
constraint_mask(::Joint{T,0}) where T = szeros(T,0,3)
constraint_mask(joint::Joint{T,1}) where T = joint.axis_mask3
constraint_mask(joint::Joint{T,2}) where T = [joint.axis_mask1; joint.axis_mask2]
constraint_mask(::Joint{T,3}) where T = SMatrix{3,3,T,9}(I)

nullspace_mask(::Joint{T,0}) where T = SMatrix{3,3,T,9}(I)
nullspace_mask(joint::Joint{T,1}) where T = [joint.axis_mask1; joint.axis_mask2]
nullspace_mask(joint::Joint{T,2}) where T = joint.axis_mask3
nullspace_mask(::Joint{T,3}) where T = szeros(T,0,3)

# impulse maps
impulse_map(relative::Symbol, joint::Joint, pbody::Node, cbody::Node, λ) = impulse_map(relative, joint, current_configuration(pbody.state)..., current_configuration(cbody.state)..., λ)

function impulse_map(relative::Symbol, joint::Joint{T,Nλ,Nb},
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion,
        η) where {T,Nλ,Nb}
    return impulse_transform(relative, joint, xa, qa, xb, qb) * impulse_projector(joint)
end

function impulse_map(relative::Symbol, joint::Joint{T,Nλ,0},
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion,
    η) where {T,Nλ}
    J = constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)
    G = cat(Diagonal(sones(3)), LVᵀmat(relative == :parent ? qa : qb), dims=(1,2))
    return Diagonal([sones(3); 0.5 * sones(3)]) * transpose(J * G)
end

function impulse_projector(joint::Joint{T,Nλ,0}) where {T,Nλ}
    zerodimstaticadjoint(constraint_mask(joint))
end

function impulse_projector(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb}
    zerodimstaticadjoint([szeros(Nb,3); -nullspace_mask(joint); nullspace_mask(joint); constraint_mask(joint)])
end

# inputs
function set_input!(joint::Joint, input::SVector)
    joint.input = zerodimstaticadjoint(nullspace_mask(joint)) * input
    return
end

set_input!(joint::Joint) = return

function add_input!(joint::Joint, input::SVector)
    joint.input += zerodimstaticadjoint(nullspace_mask(joint)) * input
    return
end

add_input!(joint::Joint) = return

function input_jacobian_control(relative::Symbol, joint::Joint, pbody::Node, cbody::Node)
    return input_jacobian_control(relative, joint, current_configuration(pbody.state)..., current_configuration(cbody.state)...) * zerodimstaticadjoint(nullspace_mask(joint))
end

# minimal coordinates
minimal_coordinates(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 3 - Nλ)

function minimal_coordinates(joint::Joint, pbody::Node, cbody::Node)
    return minimal_coordinates(joint, current_configuration(pbody.state)..., current_configuration(cbody.state)...)
end

function minimal_velocities(joint::Joint, pnode::Node, cnode::Node, timestep)
	minimal_velocities(joint, initial_configuration_velocity(pnode.state)...,
		initial_configuration_velocity(cnode.state)..., timestep)
end

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

input_dimension(joint::Joint{T,N}) where {T,N} = 3 - N
