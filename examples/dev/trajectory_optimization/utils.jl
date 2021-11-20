function step!(mech::Mechanism, x, u;
    btol=1.0e-4, undercut=Inf,
    control_inputs = (a, b, c) -> nothing)

    # set data
    data = [x; u]

    off = 0

    for body in mech.bodies
        x2, v15, q2, ω15 = unpackdata(data[off+1:end]); off += 13
        body.state.x1 = x2 - v15 * mech.Δt
        body.state.v15 = v15
        body.state.q1 = UnitQuaternion(q2...) * ωbar(-ω15, mech.Δt) * mech.Δt / 2.0
        body.state.ϕ15 = ω15
    end

    # controller
    controller!(mech, k) = control_inputs(mech, k, u)

    # simulate
    storage = simulate!(mech, mech.Δt,
        controller!, record=true, verbose=false, solver=:mehrotra!, btol=btol, undercut=undercut)

    # next state
    x_next = []

    for body in mech.bodies
        x3 = body.state.x2[1]
        v25 = body.state.vsol[2]
        _q3 = body.state.q2[1]
        q3 = [_q3.w; _q3.x; _q3.y; _q3.z]
        ω25 = body.state.ϕsol[2]
        push!(x_next, [x3; v25; q3; ω25]...)
    end

    return x_next
end

function step_grad_x!(mech::Mechanism, x, u;
    btol=1.0e-3, undercut=1.5,
    control_inputs = (a, b, c) -> nothing)
    m = length(u)
    ∂step∂x = fdjac(w -> step!(mech, w[1:(end-m)], w[(end-m+1):end],
        btol=btol, undercut=undercut, control_inputs=control_inputs),
        [x; u])[:, 1:(end-m)]

    return ∂step∂x
end

function step_grad_u!(mech::Mechanism, x, u;
    btol=1.0e-3, undercut=1.5,
    control_inputs = (a, b, c) -> nothing)
    m = length(u)
    ∂step∂u = fdjac(w -> step!(mech, w[1:(end-m)], w[(end-m+1):end],
        btol=btol, undercut=undercut, control_inputs=control_inputs),
        [x; u])[:, (end-m+1):end]

    return ∂step∂u
end

function generate_storage(mech, x)
    steps = length(x)
    nbodies = length(mech.bodies)
    storage = Storage{Float64}(steps, nbodies)

    for t = 1:steps
        off = 0
        for (i, body) in enumerate(mech.bodies)
            storage.x[i][t] = x[t][off .+ (1:3)]
            storage.v[i][t] = x[t][off .+ (4:6)]
            storage.q[i][t] = UnitQuaternion(x[t][off .+ (7:10)]...)
            storage.ω[i][t] = x[t][off .+ (11:13)]
            off += 13
        end
    end

    return storage
end





function fast_grad_x!(mech::Mechanism{T,Nn,Ne,Nb}, x, u; btol=1.0e-3, undercut=1.5,
        control_inputs = (a, b, c) -> nothing) where {T,Nn,Ne,Nb}

    step!(mech, x, u; btol = btol, undercut = undercut, control_inputs = control_inputs)

    ∂step∂x = full_matrix(mech.system) * full_data_matrix(mech)[:,1:12Nb]
    return ∂step∂x
end


function fast_grad_u!(mech::Mechanism{T,Nn,Ne,Nb}, x, u; btol=1.0e-3, undercut=1.5,
        control_inputs = (a, b, c) -> nothing) where {T,Nn,Ne,Nb}

    step!(mech, x, u; btol = btol, undercut = undercut, control_inputs = control_inputs)

    ∂step∂u = full_matrix(mech.system) * full_data_matrix(mech)[:,1:12Nb]
    return ∂step∂u
end


#
# using BenchmarkTools
# # @benchmark
# full_data_matrix(mech)
# @benchmark step_grad_x!(mech, z[1], u_control)
# @benchmark fast_grad_x!(mech, z[1], u_control)
