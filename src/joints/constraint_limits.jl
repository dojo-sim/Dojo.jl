@inline function constraint(joint::Joint{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    vertices = joint.vertices

    e1 = unlimited_constraint(joint, xa, qa, xb, qb, η)
    e2 = minimal_coordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[1:Nb½] - (joint.joint_limits[2] - e2);
            s[Nb½ .+ (1:Nb½)] - (e2 - joint.joint_limits[1]);
            constraint_mask(joint) * e1;
           ]
end

@inline function constraint_jacobian_configuration(joint::Joint{T,Nλ,Nb,N}, η) where {T,Nλ,Nb,N}
    s, γ = get_sγ(joint, η)

    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline function constraint_jacobian_parent(joint::Joint, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    # X = FiniteDiff.finite_difference_jacobian(x -> g(joint, x, qa, xb, qb, η), xa)
    # Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, UnitQuaternion(q..., false), xb, qb, η), vector(qa))
    ∇comp = szeros(T,Nb,6)
    ∇mincoord = minimal_coordinates_jacobian_configuration(:parent, joint, xa, qa, xb, qb)
    ∇unlim = unlimited_constraint_jacobian_parent(joint, xa, qa, xb, qb, η)
    return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
end

@inline function constraint_jacobian_child(joint::Joint, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    # X = FiniteDiff.finite_difference_jacobian(x -> g(joint, xa, qa, x, qb, η), xb)
    # Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, UnitQuaternion(q..., false), η), vector(qb))
    ∇comp = szeros(T,Nb,6)
    ∇mincoord = minimal_coordinates_jacobian_configuration(:child, joint, xa, qa, xb, qb)
    ∇unlim = unlimited_constraint_jacobian_child(joint, xa, qa, xb, qb, η)
    return [∇comp; ∇mincoord; -∇mincoord; ∇eqc]
end
