@inline get_body(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = collect(mechanism.bodies)[id-Ne]
@inline get_body(mechanism::Mechanism, id::Nothing) = mechanism.origin

function get_body(mechanism::Mechanism, name::Symbol)
    if name == :origin
        return mechanism.origin
    else
        for body in mechanism.bodies
            if body.name == name
                return body
            end
        end
    end
    return
end

@inline get_joint_constraint(mechanism::Mechanism, id::Integer) = mechanism.eqconstraints[id]

function get_joint_constraint(mechanism::Mechanism, name::Symbol)
    for eqc in mechanism.eqconstraints
        if eqc.name == name
            return eqc
        end
    end
    return
end

@inline get_contact_constraint(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = mechanism.ineqconstraints[id-Ne-Nb]
function get_contact_constraint(mechanism::Mechanism, name::Symbol)
    for ineqc in mechanism.ineqconstraints
        if ineqc.name == name
            return ineqc
        end
    end
    return
end

function get_node(mechanism::Mechanism{T,Nn,Ne,Nb}, id::Integer) where {T,Nn,Ne,Nb}
    if id <= Ne
        return get_joint_constraint(mechanism, id)
    elseif id <= Ne+Nb
        return get_body(mechanism, id)
    else
        return get_contact_constraint(mechanism, id)
    end
end

function get_node(mechanism::Mechanism, name::Symbol)
    node = get_body(mechanism,name)
    if node === nothing
        node = get_joint_constraint(mechanism,name)
    end
    if node === nothing
        node = get_contact_constraint(mechanism,name)
    end
    return node
end

@inline function discretize_state!(mechanism::Mechanism)
    for body in mechanism.bodies 
        discretize_state!(body, mechanism.Δt) 
    end
    return
end

@inline function ∂gab∂ʳba(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body1::Body, body2::Body) where {T,Nn,Ne,Nb,Ni}
    Δt = mechanism.Δt
    _, _, q1, ω1 = current_configuration_velocity(body1.state)
    _, _, q2, ω2 = current_configuration_velocity(body2.state)

    x1, q1 = next_configuration(body1.state, Δt)
    x2, q2 = next_configuration(body2.state, Δt)

    dGab = szeros(6,6)
    dGba = szeros(6,6)

    for connectionid in connections(mechanism.system, body1.id)
        !(connectionid <= Ne) && continue # body
        eqc = get_node(mechanism, connectionid)
        Nc = length(eqc.childids)
        off = 0
        if body1.id == eqc.parentid
            for i in 1:Nc
                joint = eqc.constraints[i]
                Nj = length(joint)
                if body2.id == eqc.childids[i]
                    Aᵀ = zerodimstaticadjoint(constraintmat(joint))
                    eqc.isspring && (dGab -= ∂springforcea∂velb(joint, body1, body2, Δt)) #should be useless
                    eqc.isdamper && (dGab -= ∂damperforcea∂velb(joint, body1, body2, Δt))
                    eqc.isspring && (dGba -= ∂springforceb∂vela(joint, body1, body2, Δt)) #should be useless
                    eqc.isdamper && (dGba -= ∂damperforceb∂vela(joint, body1, body2, Δt))
                end
                off += Nj
            end
        elseif body2.id == eqc.parentid
            for i = 1:Nc
                joint = eqc.constraints[i]
                Nj = length(joint)
                if body1.id == eqc.childids[i]
                    Aᵀ = zerodimstaticadjoint(constraintmat(joint))
                    # eqc.isspring && (dGab -= ∂springforcea∂velb(joint, body2, body1, Δt)) #should be useless
                    eqc.isdamper && (dGab -= ∂damperforcea∂velb(joint, body2, body1, Δt))
                    # eqc.isspring && (dGba -= ∂springforceb∂vela(joint, body2, body1, Δt)) #should be useless
                    eqc.isdamper && (dGba -= ∂damperforceb∂vela(joint, body2, body1, Δt))
                end
                off += Nj
            end
        end
    end
    return dGab, dGba
end
