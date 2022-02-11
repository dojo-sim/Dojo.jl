@inline function constraint(joint::Joint{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = unlimited_constraint(joint, xa, qa, xb, qb, η)
    e2 = minimal_coordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[SUnitRange(1,Nb½)] - (joint.joint_limits[2] .- e2);
            s[SUnitRange(Nb½+1,Nb)] - (e2 .- joint.joint_limits[1]);
            e1;
           ]
end

@inline function constraint_jacobian(jacobian_relative::Symbol, joint::Joint,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    (jacobian_relative == :parent) && (return constraint_jacobian_parent(joint, xa, qa, xb, qb, η))
    (jacobian_relative == :child) && (return constraint_jacobian_child(joint, xa, qa, xb, qb, η))
end

# limited
@inline function constraint_jacobian_parent(joint::Joint{T,Nλ,Nb}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb}
    ∇comp = szeros(T,Nb,7)
    ∇mincoord = minimal_coordinates_jacobian_configuration(:parent, joint, xa, qa, xb, qb, attjac=false)
    ∇unlim = unlimited_constraint_jacobian_parent(joint, xa, qa, xb, qb, η)
    return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
end
@inline function constraint_jacobian_child(joint::Joint{T,Nλ,Nb}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb}
    ∇comp = szeros(T,Nb,7)
    ∇mincoord = minimal_coordinates_jacobian_configuration(:child, joint, xa, qa, xb, qb, attjac=false)
    ∇unlim = unlimited_constraint_jacobian_child(joint, xa, qa, xb, qb, η)
    return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
end

# Unlimited
@inline function constraint_jacobian_parent(joint::Joint{T,Nλ,0}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    unlimited_constraint_jacobian_parent(joint, xa, qa, xb, qb, η)
end
@inline function constraint_jacobian_child(joint::Joint{T,Nλ,0}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    unlimited_constraint_jacobian_child(joint, xa, qa, xb, qb, η)
end
