function unpackdata(data::AbstractVector)
    x2 = data[SVector{3,Int}(1:3)]
    v15 = data[SVector{3,Int}(4:6)]
    q2 = data[SVector{4,Int}(7:10)]
    ϕ15 = data[SVector{3,Int}(11:13)]
    return x2, v15, q2, ϕ15
end

function setdata!(mechanism::Mechanism, data::AbstractVector)
    Δt = mechanism.Δt
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpackdata(data[off+1:end]); off += 13
        q2 = UnitQuaternion(q2..., false)
        body.state.x1 = getx3(x2, -v15, Δt)
        body.state.v15 = v15
        body.state.q1 = getq3(q2, -ϕ15, Δt)
        body.state.ϕ15 = ϕ15
        body.state.x2[1] = x2
        body.state.q2[1] = q2
        # setsolution!(body)
        body.state.F2[1] = SVector{3}([0,0,0.])
        body.state.τ2[1] = SVector{3}([0,0,0.])
    end
    for eqc in mechanism.eqconstraints
        dim = controldim(eqc)
        if dim > 0
            u = data[off .+ (1:dim)]; off += dim
            setForce!(eqc, u)
        end
    end

    for c in mechanism.eqconstraints 
        applyFτ!(c, mechanism, false) 
    end
    return nothing
end

function getdata(mechanism::Mechanism{T}) where T
    data = Vector{T}()
    for body in mechanism.bodies
        x2 = body.state.x2[1]
        v15 = body.state.v15
        q2 = vector(body.state.q2[1])
        ϕ15 = body.state.ϕ15
        push!(data, [x2; v15; q2; ϕ15]...)
    end
    for eqc in mechanism.eqconstraints
        if controldim(eqc) > 0
            tra = eqc.constraints[findfirst(x -> typeof(x) <: Translational, eqc.constraints)]
            rot = eqc.constraints[findfirst(x -> typeof(x) <: Rotational, eqc.constraints)]
            F = tra.Fτ
            τ = rot.Fτ
            u = [nullspacemat(tra) * F; nullspacemat(rot) * τ]
            push!(data, u...)
        end
    end
    return data
end

function setsolution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
    off = 0
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        nλ = length(eqc)
        λ = sol[off .+ (1:nλ)]; off += nλ
        eqc.λsol[2] = λ
    end
    for (i,body) in enumerate(mechanism.bodies)
        nv = 3
        nω = 3
        v25 = sol[off .+ (1:nv)]; off += nv
        ϕ25 = sol[off .+ (1:nω)]; off += nω
        body.state.vsol[2] = v25
        body.state.ϕsol[2] = ϕ25
    end
    for (i,ineqc) in enumerate(mechanism.ineqconstraints)
        N = length(ineqc)
        N½ = Int(N/2)
        s = sol[off .+ (1:N½)]; off += N½
        γ = sol[off .+ (1:N½)]; off += N½
        ineqc.ssol[2] = s
        ineqc.γsol[2] = γ
    end
    return nothing
end

function getsolution(mechanism::Mechanism{T}) where T
    sol = T[]
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        λ = eqc.λsol[2]
        push!(sol, λ...)
    end
    for (i,body) in enumerate(mechanism.bodies)
        v25 = body.state.vsol[2]
        ϕ25 = body.state.ϕsol[2]
        push!(sol, [v25; ϕ25]...)
    end
    for (i,ineqc) in enumerate(mechanism.ineqconstraints)
        s = ineqc.ssol[2]
        γ = ineqc.γsol[2]
        push!(sol, [s; γ]...)
    end
    return sol
end

function evaluate_solution!(mechanism::Mechanism, data::AbstractVector; ϵr=1e-8, ϵb=1.0e-8)
    setdata!(mechanism, data)
    mehrotra!(mechanism, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut=1.2, verbose = true))
    sol = getsolution(mechanism)
    return sol
end

function evaluate_residual!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
    system = mechanism.system
    setdata!(mechanism, data)
    setsolution!(mechanism, sol)
    setentries!(mechanism)
    return full_vector(system)
end

function finitediff_sensitivity(mechanism::Mechanism, data::AbstractVector; ϵr = 1e-8, ϵb=1.0e-8, δ = 1e-5, verbose = false)
    ndata = datadim(mechanism, attjac = false)
    nsol = soldim(mechanism)
    jac = zeros(nsol, ndata)

    for i = 1:ndata
        verbose && println("$i / $ndata")
        datap = deepcopy(data)
        datam = deepcopy(data)
        datap[i] += δ
        datam[i] -= δ

        mechanismp = deepcopy(mechanism)
        setdata!(mechanismp, deepcopy(datap))
        mehrotra!(mechanismp, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solp = getsolution(mechanismp)

        mechanismm = deepcopy(mechanism)
        setdata!(mechanismm, deepcopy(datam))
        mehrotra!(mechanismm, opts = InteriorPointOptions(rtol = ϵr, btol = ϵb, undercut = 1.2, verbose = false))
        solm = getsolution(mechanismm)

        jac[:,i] = (solp - solm) / (2δ)
    end
    return jac
end

function finitediff_data_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = soldim(mechanism)
    ndata = datadim(mechanism, attjac = false)
    jac = zeros(nsol, ndata)

    setdata!(mechanism, deepcopy(data))
    setsolution!(mechanism, deepcopy(sol))

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
    nsol = soldim(mechanism)
    jac = zeros(nsol, nsol)

    setdata!(mechanism, data)
    setsolution!(mechanism, sol)

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
