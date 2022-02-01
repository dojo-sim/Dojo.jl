function unpack_data(data::AbstractVector)
    x2 = data[SVector{3,Int}(1:3)]
    v15 = data[SVector{3,Int}(4:6)]
    q2 = data[SVector{4,Int}(7:10)]
    ϕ15 = data[SVector{3,Int}(11:13)]
    return x2, v15, q2, ϕ15
end

function set_data!(mechanism::Mechanism, data::AbstractVector)
    timestep = mechanism.timestep
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpack_data(data[off+1:end]); off += 13
        q2 = UnitQuaternion(q2..., false)
        body.state.x1 = next_position(x2, -v15, timestep)
        body.state.v15 = v15
        body.state.q1 = next_orientation(q2, -ϕ15, timestep)
        body.state.ϕ15 = ϕ15
        body.state.x2[1] = x2
        body.state.q2[1] = q2
        # set_solution!(body)
        body.state.F2[1] = SVector{3}([0,0,0.])
        body.state.τ2[1] = SVector{3}([0,0,0.])
    end
    for joint in mechanism.joints
        dim = control_dimension(joint)
        if dim > 0
            u = data[off .+ (1:dim)]; off += dim
            set_input!(joint, u)
        end
    end

    for c in mechanism.joints
        apply_input!(c, mechanism, false)
    end
    return nothing
end

function get_data(mechanism::Mechanism{T}) where T
    data = Vector{T}()
    for body in mechanism.bodies
        x2 = body.state.x2[1]
        v15 = body.state.v15
        q2 = vector(body.state.q2[1])
        ϕ15 = body.state.ϕ15
        push!(data, [x2; v15; q2; ϕ15]...)
    end
    for joint in mechanism.joints
        if control_dimension(joint) > 0
            tra = joint.constraints[findfirst(x -> typeof(x) <: Translational, joint.constraints)]
            rot = joint.constraints[findfirst(x -> typeof(x) <: Rotational, joint.constraints)]
            F = tra.Fτ
            τ = rot.Fτ
            u = [nullspace_mask(tra) * F; nullspace_mask(rot) * τ]
            push!(data, u...)
        end
    end
    return data
end

function set_solution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
    off = 0
    for (i,joint) in enumerate(mechanism.joints)
        nλ = length(joint)
        λ = sol[off .+ (1:nλ)]; off += nλ
        joint.λsol[2] = λ
    end
    for (i,body) in enumerate(mechanism.bodies)
        nv = 3
        nω = 3
        v25 = sol[off .+ (1:nv)]; off += nv
        ϕ25 = sol[off .+ (1:nω)]; off += nω
        body.state.vsol[2] = v25
        body.state.ϕsol[2] = ϕ25
    end
    for (i,contact) in enumerate(mechanism.contacts)
        N = length(contact)
        N½ = Int(N/2)
        s = sol[off .+ (1:N½)]; off += N½
        γ = sol[off .+ (1:N½)]; off += N½
        contact.ssol[2] = s
        contact.γsol[2] = γ
    end
    return nothing
end

function get_solution(mechanism::Mechanism{T}) where T
    sol = T[]
    for (i,joint) in enumerate(mechanism.joints)
        λ = joint.λsol[2]
        push!(sol, λ...)
    end
    for (i,body) in enumerate(mechanism.bodies)
        v25 = body.state.vsol[2]
        ϕ25 = body.state.ϕsol[2]
        push!(sol, [v25; ϕ25]...)
    end
    for (i,contact) in enumerate(mechanism.contacts)
        s = contact.ssol[2]
        γ = contact.γsol[2]
        push!(sol, [s; γ]...)
    end
    return sol
end

function evaluate_residual!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
    system = mechanism.system
    set_data!(mechanism, data)
    set_solution!(mechanism, sol)
    set_entries!(mechanism)
    return full_vector(system)
end

function finitediff_sensitivity(mechanism::Mechanism, data::AbstractVector; ϵr = 1e-8, ϵb=1.0e-8, δ = 1e-5, verbose = false)
    ndata = data_dimension(mechanism, attjac = false)
    nsol = solution_dimension(mechanism)
    jac = zeros(nsol, ndata)

    for i = 1:ndata
        verbose && println("$i / $ndata")
        datap = deepcopy(data)
        datam = deepcopy(data)
        datap[i] += δ
        datam[i] -= δ

        mechanismp = deepcopy(mechanism)
        set_data!(mechanismp, deepcopy(datap))
        mehrotra!(mechanismp, opts = SolverOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solp = get_solution(mechanismp)

        mechanismm = deepcopy(mechanism)
        set_data!(mechanismm, deepcopy(datam))
        mehrotra!(mechanismm, opts = SolverOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solm = get_solution(mechanismm)

        jac[:,i] = (solp - solm) / (2δ)
    end
    return jac
end

function finitediff_data_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = solution_dimension(mechanism)
    ndata = data_dimension(mechanism, attjac = false)
    jac = zeros(nsol, ndata)

    set_data!(mechanism, deepcopy(data))
    set_solution!(mechanism, deepcopy(sol))

    for i = 1:ndata
        verbose && println("$i / $ndata")
        datap = deepcopy(data)
        datam = deepcopy(data)
        datap[i] += δ
        datam[i] -= δ
        rp = evaluate_residual!(deepcopy(mechanism), datap, deepcopy(sol))
        rm = evaluate_residual!(deepcopy(mechanism), datam, deepcopy(sol))
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
end

function finitediff_sol_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = solution_dimension(mechanism)
    jac = zeros(nsol, nsol)

    set_data!(mechanism, data)
    set_solution!(mechanism, sol)

    for i = 1:nsol
        verbose && println("$i / $nsol")
        solp = deepcopy(sol)
        solm = deepcopy(sol)
        solp[i] += δ
        solm[i] -= δ
        rp = evaluate_residual!(deepcopy(mechanism), data, solp)
        rm = evaluate_residual!(deepcopy(mechanism), data, solm)
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
end
