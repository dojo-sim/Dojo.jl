# joints
function complementarity(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}

    λ = joint.impulses[2][joint_impulse_index(joint, 1)]
    s, γ = split_impulses(joint.translational, λ)
    ct = s .* γ

    λ = joint.impulses[2][joint_impulse_index(joint, 2)]
    s, γ = split_impulses(joint.rotational, λ)
    cr = s .* γ

    return [ct; cr]
end

# contacts
complementarity(mechanism, contact::ContactConstraint) = contact.impulses[2] .* contact.impulses_dual[2]

function complementarity(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:NonlinearContact{T,N}}
    γ = contact.impulses[2]
    s = contact.impulses_dual[2]
    return vcat(γ[1] * s[1], cone_product(γ[SA[2,3,4]], s[SA[2,3,4]]))
end

complementarityμ(mechanism, contact::ContactConstraint) = complementarity(mechanism, contact) - mechanism.μ * neutral_vector(contact.model)
