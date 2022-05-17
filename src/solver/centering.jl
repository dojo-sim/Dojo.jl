function centering!(mechanism::Mechanism, αaff::T) where T
    system = mechanism.system

    n = 0
    ν = 0.0
    νaff = 0.0
    for contact in mechanism.contacts
        ν, νaff, n = centering!(ν, νaff, n, mechanism, contact, get_entry(system, contact.id), αaff)
    end

    for joint in mechanism.joints
        ν, νaff, n = centering!(ν, νaff, n, mechanism, joint, get_entry(system, joint.id), αaff)
    end

    ν /= n
    νaff /= n
    return ν, νaff
end

function centering!(ν, νaff, n, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, αaff) where {T,N,Nc,Cs,N½}
    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    ν += dot(s, γ)
    νaff += dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    n += cone_degree(contact)
    return ν, νaff, n
end

function centering!(ν, νaff, n, mechanism, joint::JointConstraint{T,N,Nc}, vector_entry::Entry, αaff) where {T,N,Nc}
    for (i, element) in enumerate((joint.translational, joint.rotational))
        s, γ = split_impulses(element, joint.impulses[2][joint_impulse_index(joint,i)])
        Δs, Δγ = split_impulses(element, vector_entry.value[joint_impulse_index(joint,i)])
        ν += dot(s, γ)
        νaff += dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
        n += length(s)
    end
    return ν, νaff, n
end
