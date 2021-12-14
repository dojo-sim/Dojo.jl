function center_of_mass(mechanism::Mechanism{T}, storage::Storage{T,N}, t::Int) where {T,N}
    r = zeros(T, 3)
    for (i,body) in enumerate(mechanism.bodies)
        r += body.m * storage.x[i][t]
    end
    return r ./ total_mass(mechanism)
end

function total_mass(mechanism::Mechanism{T}) where T
    w = 0.0
    for body in mechanism.bodies
        w += body.m
    end
    return w
end

"""
    Linear and angular momentum of a body using Legendre transform.
"""
function momentum_body_new(mechanism::Mechanism{T}, body::Body{T}) where {T}
    Δt = mechanism.Δt
 
    state = body.state
    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)
    
    ezg = SA{T}[0; 0; -mechanism.g]
    D1x = - 1/Δt * body.m * (x2 - x1) + Δt/2 * body.m * ezg
    D2x =   1/Δt * body.m * (x3 - x2) + Δt/2 * body.m * ezg
    D1q =   2/Δt * LVᵀmat(q1)' * Tmat() * Rmat(q2)' * Vᵀmat() * body.J * Vmat() * Lmat(q1)' * vector(q2)
    D2q =   2/Δt * LVᵀmat(q3)' * Lmat(q2) * Vᵀmat() * body.J * Vmat() * Lmat(q2)' * vector(q3)


    p_linear_body = D2x - 0.5 * state.F2[1]
    p_angular_body = D2q - 0.5 * state.τ2[1]

    α = -1.0
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        if body.id ∈ [eqc.parentid; eqc.childids]

            λ = []
            for (i, joint) in enumerate(eqc.constraints)
                λi = eqc.λsol[2][λindex(eqc, i)]
                si, γi = get_sγ(joint, λi)
                push!(λ, λi[1:length(joint)])
                push!(λ, γi)
            end

            f_joint = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * vcat(λ...)  # computed at 1.5
            eqc.isspring && (f_joint += springforce(mechanism, eqc, body)) # computed at 1.5
            eqc.isdamper && (f_joint += damperforce(mechanism, eqc, body)) # computed at 1.5

            p_linear_body += α * 0.5 * f_joint[1:3] # TODO should't we divide F by Δt to get a force instead of a force impulse
            p_angular_body += α * 0.5 * f_joint[4:6] # TODO should't we divide F by Δt to get a force instead of a force impulse
        end
    end

    p1 = [p_linear_body; rotation_matrix(q2) * p_angular_body]
    return p1
end

function momentum(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,Ns}, t::Int) where {T,Nn,Ne,Nb,Ni,Ns}
    p = zeros(T, 6)
    com = center_of_mass(mechanism, storage, t)
    mass = total_mass(mechanism)
    p_linear_body = [] 
    p_angular_body = [] 
    for body in mechanism.bodies 
        p = momentum_body_new(mechanism, body)
        push!(p_linear_body, p[1:3])
        push!(p_angular_body, p[4:6])
    end

    p_linear = sum(p_linear_body)
    p_angular = zeros(T, 3)
    v_com = p_linear ./ mass
    for (i, body) in enumerate(mechanism.bodies)
        r = body.state.x2[1] - com
        v_body = p_linear_body[i] ./ body.m
        p_angular += p_angular_body[i]
        p_angular += cross(r, body.m * (v_body - v_com))
    end

    return [p_linear; p_angular] # in world frame
end

function momentum(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    m = [szeros(T,6) for i = 1:N]
    for i = 1:N
        m[i] = momentum(mechanism, storage, i)
    end
    return m # in world frame
end

"""
    Total mechanical energy of a mechanism.
"""
function mechanicalEnergy(mechanism::Mechanism, storage::Storage)
    ke = kineticEnergy(mechanism, storage)
    pe = potentialEnergy(mechanism, storage)
    me = ke + pe
    return me
end

"""
    Kinetic energy of a mechanism due to linear and angular velocity.
"""
function kineticEnergy(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    ke = zeros(T,N)
    for i = 1:N
        ke[i] = kineticEnergy(mechanism, storage, i)
    end
    return ke
end

"""
    Kinetic energy of a mechanism due to linear and angular velocity.
"""
function kineticEnergy(mechanism::Mechanism, storage::Storage{T,N}, t::Int) where {T,N}
    ke = 0.0
    for (i,body) in enumerate(mechanism.bodies)
        vl = storage.vl[i][t]
        ωl = storage.ωl[i][t]
        ke += 0.5 * body.m * vl' * vl
        ke += 0.5 * ωl' * body.J * ωl 
    end
    return ke
end

"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potentialEnergy(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    pe = zeros(T,N)
    for i = 1:N
        pe[i] = potentialEnergy(mechanism, storage, i)
    end
    return pe
end

"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potentialEnergy(mechanism::Mechanism{T,Nn,Ne,Nb}, storage::Storage{T,Ns}, t::Int) where {T,Nn,Ne,Nb,Ns}
    pe = 0.0
    # Gravity
    for (i,body) in enumerate(mechanism.bodies)
        x = storage.x[i][t] # x2
        q = storage.q[i][t] # q2
        z = x[3]
        pe += - body.m * mechanism.g * z
    end

    # Springs
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        if eqc.isspring
            for (j,joint) in enumerate(eqc.constraints)
                # Child
                childid = eqc.childids[j]
                xb = storage.x[childid - Ne][t] # TODO this is sketchy way to get the correct index
                qb = storage.q[childid - Ne][t] # TODO this is sketchy way to get the correct index
                # Parent
                parentid = eqc.parentid
                if parentid != nothing
                    xa = storage.x[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                    qa = storage.q[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                    (typeof(joint) <: Translational) && (force = springforceb(joint, xa, qa, xb, qb)) # actual force not impulse
                    (typeof(joint) <: Rotational) && (w = gq(joint, qa, qb).w) # actual force not impulse
                else
                    (typeof(joint) <: Translational) && (force = springforceb(joint, xb, qb)) # actual force not impulse
                    (typeof(joint) <: Rotational) && (w = gq(joint, qb).w) # actual force not impulse
                end
                # @show force
                spring = joint.spring
                if spring > 0
                    (typeof(joint) <: Translational) && (pe += 0.5 * force' * force ./ spring)
                    (typeof(joint) <: Rotational) && (pe += spring * (1 - w^2))
                end
            end
        end
    end
    return pe
end