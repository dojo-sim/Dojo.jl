"""
    getbody!(mechanism, id)

Gets the body with ID `id` from `mechanism` if it exists. If `id = nothing`, the origin will be returned.
"""
@inline getbody(mechanism::Mechanism, id::Integer) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin

"""
    getbody!(mechanism, name)

Gets the body with name `name` from `mechanism` if it exists.
"""
function getbody(mechanism::Mechanism, name::String)
    if mechanism.origin.name == name
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

"""
    geteqconstraint!(mechanism, id)

Gets the equality constraint with ID `id` from `mechanism` if it exists.
"""
@inline geteqconstraint(mechanism::Mechanism, id::Integer) = mechanism.eqconstraints[id]

"""
    geteqconstraint!(mechanism, name)

Gets the equality constraint with name `name` from `mechanism` if it exists.
"""
function geteqconstraint(mechanism::Mechanism, name::String)
    for eqc in mechanism.eqconstraints
        if eqc.name == name
            return eqc
        end
    end
    return
end

@inline getineqconstraint(mechanism::Mechanism, id::Integer) = mechanism.ineqconstraints[id]
function getineqconstraint(mechanism::Mechanism, name::String)
    for ineqc in mechanism.ineqconstraints
        if ineqc.name == name
            return ineqc
        end
    end
    return
end

"""
    getcomponent(mechanism, id)

Gets the component (body or equality constraint) with ID `id` from `mechanism` if it exists.
"""
function getcomponent(mechanism::Mechanism{T,Nn,Ne,Nb}, id::Integer) where {T,Nn,Ne,Nb}
    if id <= Ne
        return geteqconstraint(mechanism, id)
    elseif id <= Ne+Nb
        return getbody(mechanism, id)
    else
        return getineqconstraint(mechanism, id)
    end
end
getcomponent(mechanism::Mechanism, id::Nothing) = mechanism.origin

"""
    getcomponent!(mechanism, name)

Gets the component (body or equality constraint) with name `name` from `mechanism` if it exists.
"""
function getcomponent(mechanism::Mechanism, name::String)
    component = getbody(mechanism,name)
    if component === nothing
        component = geteqconstraint(mechanism,name)
    end
    if component === nothing
        component = getineqconstraint(mechanism,name)
    end
    return component
end

@inline function discretizestate!(mechanism::Mechanism)
    foreach(discretizestate!, mechanism.bodies, mechanism.Δt)
    return
end

@inline function residual_violation(mechanism::Mechanism)
    mechanism.rvio = 0

    foreach(residual_violation!, mechanism.eqconstraints, mechanism)
    foreach(residual_violation!, mechanism.bodies, mechanism)
    foreach(residual_violation!, mechanism.ineqconstraints, mechanism)

    return mechanism.rvio
end

@inline function residual_violation!(component::Component, mechanism::Mechanism)
    res = g(mechanism, component)
    mechanism.rvio = max(mechanism.rvio, norm(res, Inf))
    return nothing
end

@inline function residual_violation!(ineqc::InequalityConstraint, mechanism::Mechanism)
    res = g(mechanism, ineqc)
    mechanism.rvio = max(mechanism.rvio, norm(res, Inf))
    return nothing
end


@inline function bilinear_violation(mechanism::Mechanism)
    mechanism.bvio = 0

    foreach(bilinear_violation!, mechanism.ineqconstraints, mechanism)

    return mechanism.bvio
end

function bilinear_violation!(ineqc::InequalityConstraint, mechanism::Mechanism)
    comp = complementarity(mechanism, ineqc)
    mechanism.bvio = max(mechanism.bvio, norm(comp, Inf))
    return nothing
end


@inline function ∂gab∂ʳba(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body1::Body, body2::Body) where {T,Nn,Ne,Nb,Ni}
    Δt = mechanism.Δt
    _, _, q1, ω1 = fullargssol(body1.state)
    _, _, q2, ω2 = fullargssol(body2.state)

    x1, q1 = posargsnext(body1.state, Δt)
    x2, q2 = posargsnext(body2.state, Δt)

    dGab = zeros(6,6)
    dGba = zeros(6,6)

    for connectionid in connections(mechanism.system, body1.id)
        !(connectionid <= Ne) && continue # body
        eqc = getcomponent(mechanism, connectionid)
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
                    eqc.isspring && (dGab -= ∂springforcea∂velb(joint, body2, body1, Δt)) #should be useless
                    eqc.isdamper && (dGab -= ∂damperforcea∂velb(joint, body2, body1, Δt))
                    eqc.isspring && (dGba -= ∂springforceb∂vela(joint, body2, body1, Δt)) #should be useless
                    eqc.isdamper && (dGba -= ∂damperforceb∂vela(joint, body2, body1, Δt))
                end
                off += Nj
            end
        end
    end
    return dGab, dGba
end
