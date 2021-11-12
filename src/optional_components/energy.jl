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
    J = body.J
    x2, v2, q2, ω2 = fullargssol(body.state)
    ω2 = body.state.ωsol[2]
    p_linear_body = body.m * v2  - 0.5 * [0; 0; body.m * mechanism.g * Δt] - 0.5 * body.state.Fk[1]
    p_angular_body = (Δt * sqrt(4 / Δt^2.0 - ω2' * ω2) * body.J * ω2 - Δt * skew(ω2) * (body.J * ω2) - body.state.τk[1])

    for (i, eqc) in enumerate(mechanism.eqconstraints)
        f_joint = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2]
        p_linear_body += 0.5 * f_joint[1:3]
        p_angular_body += 0.5 * f_joint[4:6]
        if body.id == eqc.parentid
            for (i,joint) in enumerate(eqc.constraints)
                cbody = getbody(mechanism, eqc.childids[i])
                eqc.isspring && (p_linear_body += 0.5 * springforcea(joint, body.state, cbody.state, Δt)[1:3])
                eqc.isspring && (p_angular_body += 0.5 * springforcea(joint, body.state, cbody.state, Δt)[4:6])
                eqc.isdamper && (p_linear_body += 0.5 * damperforcea(joint, body.state, cbody.state, Δt)[1:3])
                eqc.isdamper && (p_angular_body += 0.5 * damperforcea(joint, body.state, cbody.state, Δt)[4:6])
            end
        end
        for (i,joint) in enumerate(eqc.constraints)
            if eqc.childids[i] == body.id
                if eqc.parentid != nothing
                    pbody = getbody(mechanism, eqc.parentid)
                    eqc.isspring && (p_linear_body += 0.5 * springforceb(joint, pbody.state, body.state, Δt)[1:3])
                    eqc.isspring && (p_angular_body += 0.5 * springforceb(joint, pbody.state, body.state, Δt)[4:6])
                    eqc.isdamper && (p_linear_body += 0.5 * damperforceb(joint, pbody.state, body.state, Δt)[1:3])
                    eqc.isdamper && (p_angular_body += 0.5 * damperforceb(joint, pbody.state, body.state, Δt)[4:6])
                else
                    eqc.isspring && (p_linear_body += 0.5 * springforceb(joint, body.state, Δt)[1:3])
                    eqc.isspring && (p_angular_body += 0.5 * springforceb(joint, body.state, Δt)[4:6])
                    eqc.isdamper && (p_linear_body += 0.5 * damperforceb(joint, body.state, Δt)[1:3])
                    eqc.isdamper && (p_angular_body += 0.5 * damperforceb(joint, body.state, Δt)[4:6])
                end
            end
        end
    end

    p1 = [p_linear_body; rotation_matrix(q2) * p_angular_body / 2] # TODO maybe this is the solution
    # p1 = [p_linear_body; rotation_matrix(q2) * p_angular_body]

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
