"""
    momentum(mechanism, storage) 

    mechanism's linear and angular momentum 

    mechanism: Mechanism 
    storage: Storage 
"""
function momentum(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    m = [szeros(T,6) for i = 1:N]
    for i = 1:N
        m[i] = momentum(mechanism, storage, i)
    end
    return m # in world frame
end

function momentum(mechanism::Mechanism{T}, body::Body{T}) where T
    timestep= mechanism.timestep
    state = body.state
    mass = body.mass 
    inertia = body.inertia
    # x1, q1 = previous_configuration(state)
    x2, q2 = current_configuration(state)
    x3, q3 = next_configuration(state, timestep)

    v15 = body.state.vsol[2] # v1.5
    ω15 = body.state.ωsol[2] # ω1.5

    D2x = 1 / timestep * mass * (x3 - x2) - 0.5 * timestep * mass * mechanism.gravity
    D2q = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)
    p_linear_body = D2x - 0.5 * state.F2
    p_angular_body = D2q - 0.5 * state.τ2

    for joint in mechanism.joints
        if body.id ∈ [joint.parent_id, joint.child_id]
            f_joint = impulse_map(mechanism, joint, body) * joint.impulses[2] # computed at t = 1.5

            joint.spring && (f_joint += spring_impulses(mechanism, joint, body)) # computed at t = 1.5
            joint.damper && (f_joint += damper_impulses(mechanism, joint, body)) # computed at t = 1.5

            p_linear_body -= 0.5 * f_joint[1:3]
            p_angular_body -= 0.5 * f_joint[4:6]
        end
    end

    return [p_linear_body; rotation_matrix(q2) * p_angular_body]
end

function momentum(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,Ns}, t::Int) where {T,Nn,Ne,Nb,Ni,Ns}
    p = zeros(T, 6)
    com = center_of_mass(mechanism, storage, t)
    mass = total_mass(mechanism)
    p_linear_body = [storage.px[i][t] for i = 1:Nb] # in world frame
    p_angular_body = [storage.pq[i][t] for i = 1:Nb] # in world frame

    p_linear = sum(p_linear_body)
    p_angular = zeros(T, 3)
    v_com = p_linear ./ mass
    for (i, body) in enumerate(mechanism.bodies)
        r = storage.x[i][t] - com
        v_body = p_linear_body[i] ./ body.mass
        p_angular += p_angular_body[i]
        p_angular += cross(r, body.mass * (v_body - v_com)) #TODO maybe there is cleaner way to handle the factor 2
    end

    return [p_linear; p_angular] # in world frame
end

function center_of_mass(mechanism::Mechanism{T}, storage::Storage{T,N}, t::Int) where {T,N}
    r = zeros(T, 3)
    for (i,body) in enumerate(mechanism.bodies)
        r += body.mass * storage.x[i][t]
    end
    return r ./ total_mass(mechanism)
end

function total_mass(mechanism::Mechanism{T}) where T
    w = 0.0
    for body in mechanism.bodies
        w += body.mass
    end
    return w
end
