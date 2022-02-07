spring_parent(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_parent(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
spring_child(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_child(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
damper_parent(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_parent(joint, current_configuration(bodya.state)..., bodya.state.vsol[2],
    bodya.state.ϕsol[2], current_configuration(bodyb.state)..., bodyb.state.vsol[2], bodyb.state.ϕsol[2], unitary=unitary)
damper_child(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_child(joint, current_configuration(bodya.state)..., bodya.state.vsol[2],
    bodya.state.ϕsol[2], current_configuration(bodyb.state)..., bodyb.state.vsol[2], bodyb.state.ϕsol[2], unitary=unitary)

spring_parent(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
spring_child(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_parent(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_child(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)

spring_parent(joint::Translational3{T}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)
spring_child(joint::Translational3{T}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Spring Force
################################################################################
function spring_force(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
            xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δmincoord = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb) # in the a frame
    Fτ = spring * Aᵀ * Δmincoord # in the a frame
    return Fτ
end

function spring_force_jacobian_configuration(jacobian_relative::Symbol,
        joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇Fτ = spring * Aᵀ * - minimal_coordinates_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb)
    return ∇Fτ
end

@inline function spring_relative(relative::Symbol, joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    Fτ = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
    Fτ = impulse_transform(relative, joint, xa, qa, xb, qb) * Fτ
    return Fτ
end

@inline function spring_parent(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    # @warn "remove this"
    # error()
    spring_relative(:parent, joint, xa, qa, xb, qb; unitary=unitary)
end

@inline function spring_child(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    # @warn "remove this"
    # error()
    spring_relative(:child, joint, xa, qa, xb, qb; unitary=unitary)
end


@inline function minimal_velocities(joint::Translational, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
    Δv = vrotate(Δvw, inv(qa)) # in the a frame
    return nullspace_mask(joint) * Δv
end

function damper_force(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector,
        ωb::AbstractVector; unitary::Bool=false) where T
    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Fτ = damper * Aᵀ * -minimal_velocities(joint, xa, qa, va, ωa, xb, qb, vb, ωb) # in the a frame
    return Fτ
end

@inline function damper_relative(relative::Symbol, joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector,
        ωb::AbstractVector; unitary::Bool=false) where T
    Fτ = damper_force(joint, xa, qa, va, ωa, xb, qb, vb, ωb) # in the a frame
    Fτa = impulse_transform(relative, joint, xa, qa, xb, qb) * Fτ
    return Fτa
end
@inline function damper_parent(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector,
        ωb::AbstractVector; unitary::Bool=false) where T
    # Fτ = damper_force(joint, xa, qa, va, ωa, xb, qb, vb, ωb) # in the a frame
    # Fτa = impulse_transform_parent(joint, xa, qa, xb, qb) * Fτ
    # return Fτa
    return damper_relative(:parent, joint, xa, qa, va, ωa, xb, qb, vb, ωb)
end
@inline function damper_child(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector,
        ωb::AbstractVector; unitary::Bool=false) where T
    # Fτ = damper_force(joint, xa, qa, va, ωa, xb, qb, vb, ωb) # in the a frame
    # Fτb = impulse_transform_child(joint, xa, qa, xb, qb) * Fτ
    # return Fτb
    return damper_relative(:child, joint, xa, qa, va, ωa, xb, qb, vb, ωb)
end

spring_parent_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_parent_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuraion_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)

spring_parent_jacobian_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_parent_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_configuration_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_configuration_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_configuration_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)


@inline function spring_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector,
        qb::UnitQuaternion, timestep::T; unitary::Bool=false) where T
    Fτ = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        spring_force_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb, unitary=unitary)
    ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, Fτ)
    return timestep * ∇xq
end
function spring_parent_jacobian_configuration_parent(joint::Translational, body1::Node,
        body2::Node, timestep::T; attjac::Bool = true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    attjac && (return spring_jacobian_configuration(:parent, :parent, joint, xa, qa, xb, qb, timestep; unitary=false))

    X = FiniteDiff.finite_difference_jacobian(xa -> spring_parent(joint, xa, qa, xb, qb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_parent(joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end
function spring_parent_jacobian_configuration_child(joint::Translational, body1::Node,
        body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    attjac && (return spring_jacobian_configuration(:parent, :child, joint, xa, qa, xb, qb, timestep; unitary=false))

    X = FiniteDiff.finite_difference_jacobian(xb -> spring_parent(joint, xa, qa, xb, qb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_parent(joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end
function spring_child_jacobian_configuration_child(joint::Translational, body1::Node,
        body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    attjac && (return spring_jacobian_configuration(:child, :child, joint, xa, qa, xb, qb, timestep; unitary=false))

    X = FiniteDiff.finite_difference_jacobian(xb -> spring_child(joint, xa, qa, xb, qb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_child(joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end
function spring_child_jacobian_configuraion_parent(joint::Translational, body1::Node,
        body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    attjac && (return spring_jacobian_configuration(:child, :parent, joint, xa, qa, xb, qb, timestep; unitary=false))

    X = FiniteDiff.finite_difference_jacobian(xa -> spring_child(joint, xa, qa, xb, qb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_child(joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end

function spring_jacobian_velocity(relative::Symbol, jacobian_relative::Symbol, joint::Translational,
        body1::Node, body2::Node, timestep::T) where T
    return timestep * szeros(T, 6, 6)
end
function spring_parent_jacobian_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return spring_jacobian_velocity(:parent, :parent, joint, body1, body2, timestep)
end
function spring_parent_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return spring_jacobian_velocity(:parent, :child, joint, body1, body2, timestep)
end
function spring_child_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return spring_jacobian_velocity(:child, :child, joint, body1, body2, timestep)
end
function spring_child_configuration_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return spring_jacobian_velocity(:child, :parent, joint, body1, body2, timestep)
end


@inline function damper_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T
    Fτ = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, unitary=unitary)
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_force_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, unitary=unitary)
    ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, Fτ)
    return timestep * ∇xq
end
function damper_parent_jacobian_configuration_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_configuration(:parent, :parent, joint, xa, qa, va, ωa, xb, qb, vb, ωb, timestep; unitary=false))


    X = FiniteDiff.finite_difference_jacobian(xa -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> damper_parent(joint, xa, UnitQuaternion(qa..., false), va, ωa, xb, qb, vb, ωb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end
function damper_parent_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_configuration(:parent, :child, joint, xa, qa, va, ωa, xb, qb, vb, ωb, timestep; unitary=false))


    X = FiniteDiff.finite_difference_jacobian(xb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> damper_parent(joint, xa, qa, va, ωa, xb, UnitQuaternion(qb..., false), vb, ωb), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end
function damper_child_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_configuration(:child, :child, joint, xa, qa, va, ωa, xb, qb, vb, ωb, timestep; unitary=false))


    X = FiniteDiff.finite_difference_jacobian(xb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> damper_child(joint, xa, qa, va, ωa, xb, UnitQuaternion(qb..., false), vb, ωb), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end
function damper_child_jacobian_configuration_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_configuration(:child, :parent, joint, xa, qa, va, ωa, xb, qb, vb, ωb, timestep; unitary=false))


    X = FiniteDiff.finite_difference_jacobian(xa -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> damper_child(joint, xa, UnitQuaternion(qa..., false), va, ωa, xb, qb, vb, ωb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end

@inline function damper_jacobian_velocity(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_force_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, unitary=unitary)
    return timestep * ∇xq
end
function damper_parent_jacobian_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T, attjac::Bool=true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_velocity(:parent, :parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    V = FiniteDiff.finite_difference_jacobian(va -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), va)
    Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωa)
    return timestep * [V Ω]
end
function damper_parent_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T, attjac::Bool=true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_velocity(:parent, :child, joint, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    V = FiniteDiff.finite_difference_jacobian(vb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), vb)
    Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωb)
    return timestep * [V Ω]
end
function damper_child_configuration_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T, attjac::Bool=true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_velocity(:child, :child, joint, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    V = FiniteDiff.finite_difference_jacobian(vb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), vb)
    Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωb)
    return timestep * [V Ω]
end
function damper_child_configuration_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T, attjac::Bool=true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]
    attjac && (return damper_jacobian_velocity(:child, :parent, joint, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    V = FiniteDiff.finite_difference_jacobian(va -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), va)
    Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωa)
    return timestep * [V Ω]
end
