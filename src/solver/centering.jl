function centering!(mechanism::Mechanism, αaff::T) where T
    system = mechanism.system

    parameters = [SA{T}[0;0;0]] # ν, νaff, n
    for contact in mechanism.contacts
        centering!(parameters, contact, get_entry(system, contact.id), αaff)
    end

    for joint in mechanism.joints
        centering!(parameters, joint, get_entry(system, joint.id), αaff)
    end

    ν = parameters[1][1]/parameters[1][3]
    νaff = parameters[1][2]/parameters[1][3]
    return ν, νaff
end

function centering!(parameters, contact::ContactConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, αaff) where {T,N,Nc,Cs,N½}
    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[SUnitRange(1,N½)]
    Δγ = vector_entry.value[SUnitRange(N½+1,N)]
    ν = dot(s, γ)
    νaff = dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    n = cone_degree(contact)
    parameters[1] += SA{T}[ν;νaff;n]

    return
end

function centering!(parameters, joint::JointConstraint{T,N,Nc}, vector_entry::Entry, αaff) where {T,N,Nc}
    joint_t = joint.translational
    s, γ = split_impulses(joint_t, joint.impulses[2][joint_impulse_index(joint,1)])
    Δs, Δγ = split_impulses(joint_t, vector_entry.value[joint_impulse_index(joint,1)])
    ν = dot(s, γ)
    νaff = dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    n = length(s)
    parameters[1] += SA{T}[ν;νaff;n]

    joint_r = joint.rotational
    s, γ = split_impulses(joint_r, joint.impulses[2][joint_impulse_index(joint,2)])
    Δs, Δγ = split_impulses(joint_r, vector_entry.value[joint_impulse_index(joint,2)])
    ν = dot(s, γ)
    νaff = dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    n = length(s)
    parameters[1] += SA{T}[ν;νaff;n]

    return
end