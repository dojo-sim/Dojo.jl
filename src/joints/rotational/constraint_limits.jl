@inline function constraint(joint::Rotational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = Vmat(qa \ qb / joint.qoffset)
    e2 = minimal_coordinates(joint, qa, qb)
    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[1:Nb½] - (joint.joint_limits[2] .- e2);
            s[Nb½ .+ (1:Nb½)] - (e2 .- joint.joint_limits[1]);
            constraint_mask(joint) * e1;
           ]
end

@inline function constraint_jacobian_configuration(joint::Rotational{T,Nλ,Nb,N}, η) where {T,Nλ,Nb,N}
    s, γ = get_sγ(joint, η)

    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline function constraint_jacobian_parent(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> constraint(joint, xa, UnitQuaternion(q..., false), xb, qb, η), vector(qa))
    return [X Q]
end

@inline function constraint_jacobian_child(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> constraint(joint, xa, qa, xb, UnitQuaternion(q..., false), η), vector(qb))
    return [X Q]
end
