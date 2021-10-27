function fdjac(f, x; δ = 1e-5)
    x = Vector(x)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end

function finitediff_vel(joint::AbstractJoint, pbody::AbstractBody, cbody::AbstractBody, Δt::T,
        evalf, jacf; diff_body::Symbol = :child) where {T}

    jac0 = jacf(joint, pbody, cbody, Δt)
    function f(x)
        v2 = x[1:3]
        ω2 = x[4:6]
        if diff_body == :parent
            cstate = deepcopy(cbody.state)
            pstate = deepcopy(pbody.state)
            pstate.vsol[2] = v2
            pstate.ωsol[2] = ω2
        elseif diff_body == :child
            cstate = deepcopy(cbody.state)
            cstate.vsol[2] = v2
            cstate.ωsol[2] = ω2
            if typeof(pbody) <: Origin
                return evalf(joint, cstate, Δt)
            else
                pstate = deepcopy(pbody.state)
            end
        end
        return evalf(joint, pstate, cstate, Δt)
    end

    if diff_body == :child
        v2 = cbody.state.vsol[2]
        ω2 = cbody.state.ωsol[2]
        x = [v2; ω2]
    elseif diff_body == :parent
        v2 = pbody.state.vsol[2]
        ω2 = pbody.state.ωsol[2]
        x = [v2; ω2]
    else
        error("invalid diff_body")
    end
    jac1 = fdjac(f, x)
    return jac0, jac1
end

function finitediff_pos(joint::AbstractJoint, pbody::AbstractBody, cbody::AbstractBody, Δt::T,
        evalf, jacf; diff_body::Symbol = :child) where {T}

    jac0 = jacf(joint, pbody, cbody, Δt)
    function f(x)
        x2 = x[1:3]
        q2 = x[4:7]
        if diff_body == :parent
            cstate = deepcopy(cbody.state)
            pstate = deepcopy(pbody.state)
            pstate.xk[1] = x2
            pstate.qk[1] = UnitQuaternion(q2...)
        elseif diff_body == :child
            cstate = deepcopy(cbody.state)
            cstate.xk[1] = x2
            cstate.qk[1] = UnitQuaternion(q2...)
            if typeof(pbody) <: Origin
                return evalf(joint, cstate, Δt)
            else
                pstate = deepcopy(pbody.state)
            end
        end
        return evalf(joint, pstate, cstate, Δt)
    end

    if diff_body == :child
        x2 = cbody.state.xk[1]
        q2 = cbody.state.qk[1]
        x = [x2; [q2.w, q2.x, q2.y, q2.z]]
    elseif diff_body == :parent
        x2 = pbody.state.xk[1]
        q2 = pbody.state.qk[1]
        x = [x2; [q2.w, q2.x, q2.y, q2.z]]
    else
        error("invalid diff_body")
    end
    M = [I zeros(3,3); zeros(4,3) LVᵀmat(q2)]
    jac1 = fdjac(f, x) * M
    return jac0, jac1
end
