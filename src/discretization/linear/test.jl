@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.x1 = x1
    state.q1 = q1
    state.v15 = v1
    state.ϕ15 = ω1
    state.vsol[2] = v2
    state.ϕsol[2] = ω2

    state.x2[1] = x1 + v1*Δt
    state.q2[1] = q1 * ωbar(ω1,Δt) * Δt / 2

    return
end
