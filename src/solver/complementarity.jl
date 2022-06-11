# joints
function complementarity(mechanism, joint::JointConstraint{T,N,Nc};
    scaling::Bool=false) where {T,N,Nc}

    c = []
    for (i, element) in enumerate((joint.translational, joint.rotational))
        λi = joint.impulses[2][joint_impulse_index(joint, i)]
        si, γi = split_impulses(element, λi)
        push!(c, si .* γi)
    end

    return vcat(c...)
end

# contacts
complementarity(mechanism, contact::ContactConstraint; scaling=false) = contact.impulses[2] .* contact.impulses_dual[2]

function complementarity(mechanism, contact::RigidContactConstraint{T,N,Nc,Cs,N½};
    scaling::Bool = false) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
    γ = contact.impulses[2]
    s = contact.impulses_dual[2]
    return vcat(γ[1] * s[1], cone_product(γ[@SVector [2,3,4]], s[@SVector [2,3,4]]))
end

complementarityμ(mechanism, contact::ContactConstraint; scaling=false) = complementarity(mechanism, contact, scaling=scaling) - mechanism.μ * neutral_vector(contact.model)
