@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.xc = x1
    state.qc = q1
    state.vc = v1
    state.ωc = ω1
    state.vsol[2] = v2
    state.ωsol[2] = ω2

    state.xk[1] = x1 + v1*Δt
    state.qk[1] = q1 * ωbar(ω1,Δt) * Δt / 2

    state.xsol[2] = state.xk[1]
    state.qsol[2] = state.qk[1]

    return
end
