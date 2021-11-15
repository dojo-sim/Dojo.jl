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
    J = body.J
    x1 = body.state.xk[1]
    q1 = body.state.qk[1]
    v05 = body.state.vc # v0.5
    ω05 = body.state.ωc # ω0.5
    v15 = body.state.vsol[2] # v1.5
    ω15 = body.state.ωsol[2] # ω1.5

    p_linear_body = body.m * v15  - 0.5 * [0; 0; body.m * mechanism.g * Δt] - 0.5 * body.state.Fk[1] #BEST TODO should't we divide F by Δt to get a force instead of a force impulse
    p_angular_body = Δt * sqrt(4 / Δt^2.0 - ω15' * ω15) * body.J * ω15 + Δt * skew(ω15) * (body.J * ω15) - body.state.τk[1] # TODO should't we divide tau by Δt to get a torque instead of a torque impulse

    α = -1.0
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        if body.id ∈ [eqc.parentid; eqc.childids]
            f_joint = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2] # computed at 1.5
            eqc.isspring && (f_joint += springforce(mechanism, eqc, body)) # computed at 1.5
            eqc.isdamper && (f_joint += damperforce(mechanism, eqc, body)) # computed at 1.5

            # @show f_joint[1:3]
            p_linear_body += α * 0.5 * f_joint[1:3] # TODO should't we divide F by Δt to get a force instead of a force impulse
            p_angular_body += α * 0.5 * f_joint[4:6] # TODO should't we divide F by Δt to get a force instead of a force impulse
        end
    end

    p1 = [p_linear_body; rotation_matrix(q1) * p_angular_body/sqrt(2)] # everything in the world frame, TODO maybe this is the solution
    # p1 = [p_linear_body; rotation_matrix(q1) * p_angular_body] # everything in the world frame, TODO maybe this is the solution
    return p1
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
        v_body = p_linear_body[i] ./ body.m
        p_angular += p_angular_body[i]
        # p_angular += cross(r, body.m * (v_body - v_com))
        p_angular += cross(r, body.m * (v_body - v_com))/sqrt(2) #TODO maybe there is cleaner way to handle the factor 2
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
        x = storage.x[i][t] # x1
        q = storage.q[i][t] # q1
        z = x[3]
        pe += - body.m * mechanism.g * z
    end
    # Springs
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        if eqc.isspring
            for (j,joint) in enumerate(eqc.constraints)
                @show joint
                # Child
                childid = eqc.childids[j]
                xb = storage.x[childid - Ne][t] # TODO this is sketchy way to get the correct index
                qb = storage.q[childid - Ne][t] # TODO this is sketchy way to get the correct index
                # Parent
                parentid = eqc.parentid
                if parentid != nothing
                    xa = storage.x[parentid][t]
                    qa = storage.q[parentid][t]
                    (typeof(joint) <: Translational) && (force = springforceb(joint, xa, qa, xb, qb)) # actual force not impulse
                    (typeof(joint) <: Rotational) && (force = springforceb(joint, qa, qb)) # actual force not impulse
                else
                    (typeof(joint) <: Translational) && (force = springforceb(joint, xb, qb)) # actual force not impulse
                    (typeof(joint) <: Rotational) && (force = springforceb(joint, qb)) # actual force not impulse
                end
                # @show force[3]
                # @show body.state.xk[1][3] * joint.spring
                pe += 0.5 * force' * force ./ joint.spring
            end
        end
    end
    return pe
end

function potentialEnergy(mechanism::Mechanism)
    V = 0.0
    for body in mechanism.bodies
        V += potentialEnergy(mechanism, body)
    end
    for eqc in mechanism.eqconstraints
        V += potentialEnergy(mechanism, eqc)
    end
    return V
end


"""
    Potential energy of a joint due to the spring.
"""
function potentialEnergy(mechanism::Mechanism, eqc::EqualityConstraint)
    V = 0.0
    parentid = eqc.parentid
    Δt = mechanism.Δt
    if parentid != nothing
        pbody = getbody(mechanism, parentid)
    else
        pbody = mechanism.origin
    end
    for (i,joint) in enumerate(eqc.constraints)
        cbody = getbody(mechanism, eqc.childids[i])
        V += potentialEnergy(joint, pbody, cbody, Δt)
    end
    return V
end

"""
    Potential energy of a joint due to the spring.
"""
potentialEnergy(joint::Rotational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Rotational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0

function potentialEnergy(joint::Translational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Translational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end



















################################################################################
# OLD STUFF
################################################################################


"""
    Total linear and angular momentum of a mechanism.
"""
function momentum(mechanism::Mechanism{T}) where {T}
    p = zeros(T, 6)
    com = center_of_mass(mechanism)
    mass = total_mass(mechanism)
    p_body = Vector{T}[]

    for body in mechanism.bodies
        push!(p_body, momentum_body(mechanism, body))
    end

    p_linear = sum([p[1:3] for p in p_body])
    p_angular = zeros(T, 3)
    v_com = p_linear ./ mass
    for (i, body) in enumerate(mechanism.bodies)
        r = body.state.xk[1] - com
        v_body = p_body[i][1:3] ./ body.m
        p_angular += p_body[i][4:6]
        # p_angular += cross(r, body.m * (v_body - v_com))
        p_angular += cross(r, body.m * (v_body - v_com))/2 #TODO maybe there is cleaner way to handle the factor 2
    end

    return [p_linear; p_angular]
end

"""
    center of mass of a mechanism.
"""
function center_of_mass(mechanism::Mechanism{T}) where T
    r = zeros(T, 3)
    for body in mechanism.bodies
        r += body.m * body.state.xk[1]
    end
    return r ./ total_mass(mechanism)
end

"""
    Linear and angular momentum of a body using Legendre transform.
"""
function momentum_body(mechanism::Mechanism{T}, body::Body{T}) where {T}
    Δt = mechanism.Δt
    J = body.J
    x1 = body.state.xk[1]
    q1 = body.state.qk[1]
    v05 = body.state.vc
    ω05 = body.state.ωc
    v15 = body.state.vsol[2]
    ω15 = body.state.ωsol[2]

    p_linear_body = body.m * v15  - 0.5 * [0; 0; body.m * mechanism.g * Δt] - 0.5 * body.state.Fk[1]
    # p_angular_body = Δt * sqrt(4 / Δt^2.0 - ω05' * ω05) * body.J * ω05 - Δt * skew(ω05) * (body.J * ω05) + body.state.τk[1]
    # @show body.state.Fk[1]
    # @show body.state.τk[1]
    # p_angular_body = Δt * sqrt(4 / Δt^2.0 - ω15' * ω15) * body.J * ω15 + Δt * skew(ω15) * (body.J * ω15) - 0.5*body.state.τk[1]
    p_angular_body = Δt * sqrt(4 / Δt^2.0 - ω15' * ω15) * body.J * ω15 - Δt * skew(ω15) * (body.J * ω15) - body.state.τk[1]

    α = -1.0
    for (i, eqc) in enumerate(mechanism.eqconstraints)

        f_joint = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2] # computed at 1.5
        eqc.isspring && (f_joint += springforce(mechanism, eqc, body)) # computed at 1.5
        eqc.isdamper && (f_joint += damperforce(mechanism, eqc, body)) # computed at 1.5

        p_linear_body += α * 0.5 * f_joint[1:3]
        p_angular_body += α * 0.5 * f_joint[4:6]
    end

    p1 = [p_linear_body; rotation_matrix(q1) * p_angular_body/2] # everything in the world frame, TODO maybe this is the solution
    # p1 = [p_linear_body; rotation_matrix(q2) * p_angular_body]
    return p1
end


"""
    Kinetic energy of a mechanism due to translation and rotation velocities.
"""
function kineticEnergy(mechanism::Mechanism; linear::Bool = true, angular::Bool = true)
    Δt = mechanism.Δt
    ET = 0.0
    com = center_of_mass(mechanism)
    for body in mechanism.bodies
        x2, v2, q2, ω2 = fullargssol(body.state)
        # M = [body.m * I szeros(3,3); szeros(3,3) body.J]
        p1 = momentum_body(mechanism, body) # in world frame
        p_linear_body = rotation_matrix(inv(q2)) * p1[1:3] # in body frame             , because we need to multiply by m which in world frame
        p_angular_body = rotation_matrix(inv(q2)) * p1[4:6] # in body frame, because we need to multiply by J which in body frame
        # @show scn.(p_linear_body)
        # @show scn.(p_angular_body)
        # ET += 0.5 * p_linear_body' * p_linear_body / body.m
        # ET += 0.5 * p_angular_body' * (body.J \ p_angular_body)
        v1 = p_linear_body ./ body.m
        # @show scn.(v1)
        ω1 = body.J \ p_angular_body
        # @show scn.(ω1)
        linear && (ET += 0.5 * body.m * v1' * v1)
        # @show norm(0.5 * body.m * v1' * v1 - 0.5 * p_linear_body' * p_linear_body / body.m)
        # @show norm(0.5 * ω1' * body.J * ω1 - 0.5 * p_angular_body' * (body.J \ p_angular_body))
        # @show 0.5 * body.m * v1' * v1
        # @show 0.5 * ω1' * body.J * ω1
        angular && (ET += 0.5 * ω1' * body.J * ω1)
        # @show scn.(ω1, digits = 5)
    end
    return ET
end


"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potentialEnergy(mechanism::Mechanism)
    V = 0.0
    for body in mechanism.bodies
        V += potentialEnergy(mechanism, body)
    end
    for eqc in mechanism.eqconstraints
        V += potentialEnergy(mechanism, eqc)
    end
    return V
end

"""
    Potential energy of one body due to gravity.
"""
function potentialEnergy(mechanism::Mechanism, body::Body)
    x, q = posargsk(body.state)
    z = x[3]
    V = - body.m * mechanism.g * z
    return V
end

"""
    Potential energy of a joint due to the spring.
"""
function potentialEnergy(mechanism::Mechanism, eqc::EqualityConstraint)
    V = 0.0
    parentid = eqc.parentid
    Δt = mechanism.Δt
    if parentid != nothing
        pbody = getbody(mechanism, parentid)
    else
        pbody = mechanism.origin
    end
    for (i,joint) in enumerate(eqc.constraints)
        cbody = getbody(mechanism, eqc.childids[i])
        V += potentialEnergy(joint, pbody, cbody, Δt)
    end
    return V
end

"""
    Potential energy of a joint due to the spring.
"""
potentialEnergy(joint::Rotational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Rotational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0

function potentialEnergy(joint::Translational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Translational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

"""
    Total mechanical energy of a mechanism.
"""
function mechanicalEnergy(mechanism::Mechanism)
    ke = kineticEnergy(mechanism)
    pe = potentialEnergy(mechanism)
    e = ke + pe
    return e
end
