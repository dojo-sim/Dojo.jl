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
        p_angular += p_body[i][4:6]
        p_angular += cross(r, body.m * (p_body[i][1:3] ./ body.m - 1.0 * v_com))
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
function momentum_body(mechanism::Mechanism{T}, body::Body{T}) where {T}
    Δt = mechanism.Δt
    x2, v2, q2, ω2 = fullargssol(body.state)

    p_linear_body = body.m * v2  - 0.5 * [0; 0; body.m * mechanism.g * Δt] - 0.5 * body.state.Fk[1]
    p_angular_body = rotation_matrix(q2) * (Δt * skewplusdiag(ω2, sqrt(4 / Δt^2 - ω2' * ω2)) * (body.J * ω2)) - body.state.τk[1]

    p1 = [p_linear_body; p_angular_body]
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        # f_joint = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2]
        # p_linear_body -= 0.5 * f_joint[1:3]
        # p_angular_body -= 0.5 * f_joint[4:6]

        if body.id == eqc.parentid
            for (i,joint) in enumerate(eqc.constraints)
                cbody = getbody(mechanism, eqc.childids[i])
                # eqc.isspring && (p1 += 0.5 * springforcea(joint, body.state, cbody.state, Δt))
                # eqc.isdamper && (p1 += 0.5 * damperforcea(joint, body.state, cbody.state, Δt))
            end
        end
        for (i,joint) in enumerate(eqc.constraints)
            if eqc.childids[i] == body.id
                if eqc.parentid != nothing
                    pbody = getbody(mechanism, eqc.parentid)
                    # eqc.isspring && (p1 += 0.5 * springforceb(joint, pbody.state, body.state, Δt))
                    # eqc.isdamper && (p1 += 0.5 * damperforceb(joint, pbody.state, body.state, Δt))
                else
                    # eqc.isspring && (p1 += 0.5 * springforceb(joint, body.state, Δt))
                    # eqc.isdamper && (p1 += 0.5 * damperforceb(joint, body.state, Δt))
                end
            end
        end
    end

    # return [p_linear_body; p_angular_body]
    return p1
end


"""
    Total mechanical energy of a mechanism.
"""
function mechanicalEnergy(mechanism::Mechanism)
    ET = kineticEnergy(mechanism)
    EV = potentialEnergy(mechanism)
    E = ET + EV
    return E
end

"""
    Kinetic energy of a mechanism due to translation and rotation velocities.
"""
function kineticEnergy(mechanism::Mechanism)
    Δt = mechanism.Δt
    ET = 0.0
    com = center_of_mass(mechanism)
    for body in mechanism.bodies
        p1 = momentum_body(mechanism, body)
        M = [body.m * I szeros(3,3); szeros(3,3) body.J]
        v1 = M \ p1
        ET += 0.5 * v1'* M * v1
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
