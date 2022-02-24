function set_maximal_coordinates!(body::Body; 
    x::AbstractVector=SA[0;0;0], q::UnitQuaternion=one(UnitQuaternion))

    body.state.x2[1] = x
    body.state.q2[1] = q

    return body.state.x2[1], body.state.q2[1]
end

function set_maximal_coordinates!(pbody::Node, cbody::Body;
        parent_vertex::AbstractVector=SA[0;0;0], child_vertex::AbstractVector=SA[0;0;0],
        Δx::AbstractVector=SA[0;0;0], Δq::UnitQuaternion=one(UnitQuaternion))

    q1 = pbody.state.q2[1]
    q2 = pbody.state.q2[1] * Δq
    x2 = pbody.state.x2[1] + vector_rotate(parent_vertex + Δx, q1) - vector_rotate(child_vertex, q2)

    return set_maximal_coordinates!(cbody; x = x2, q = q2)
end

function set_maximal_velocities!(body::Body; 
    v::AbstractVector=SA[0;0;0], ω::AbstractVector=SA[0;0;0])

    body.state.v15 = v
    body.state.ϕ15 = ω

    return body.state.v15, body.state.ϕ15
end

function set_maximal_velocities!(pbody::Node, cbody::Body;
        parent_vertex::AbstractVector=SA[0;0;0], child_vertex::AbstractVector=SA[0;0;0],
        Δv::AbstractVector=SA[0;0;0], Δω::AbstractVector=SA[0;0;0])

    x1 = pbody.state.x2[1]
    v1 = pbody.state.v15
    q1 = pbody.state.q2[1]
    ω1 = pbody.state.ϕ15 # in local coordinates

    x2 = cbody.state.x2[1]
    # v2 = cbody.state.v15
    q2 = cbody.state.q2[1]
    # ω2 = cbody.state.ϕ15 # in local coordinates

    # Ω(B/W)b = Ra->b * [Ω(B/A)a + Ω(A/W)a]
    ω2 = vector_rotate(Δω + ω1, inv(q2) * q1)
    # V(cb,B/W)w =
    ω1w = vector_rotate(ω1, q1)
    ω2w = vector_rotate(ω2, q2)
    Δvw = vector_rotate(Δv, q1)
    cApB_w = (x2 + vector_rotate(child_vertex, q2)) - x1
    pBcB_w = - vector_rotate(child_vertex, q2)
    v2 = copy(v1)
    v2 += skew(ω1w) * cApB_w
    v2 += skew(ω2w) * pBcB_w
    v2 += Δvw

    return set_maximal_velocities!(cbody; v = v2, ω = ω2)
end

function set_input!(body::Body;
        F::AbstractVector=SA[0;0;0], τ::AbstractVector=SA[0;0;0], p::AbstractVector=SA[0;0;0])
    # F and p in local coordinates
    τ += torque_from_force(F, p) # in local coordinates
    set_input!(body.state, vector_rotate(F,body.state.q2[1]), τ)
    return
end

torque_from_force(F::AbstractVector, p::AbstractVector) = cross(p, F)
