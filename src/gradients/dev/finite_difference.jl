function get_solution(mechanism::Mechanism{T}) where T
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

function set_solution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
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

function evaluate_residual!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
    system = mechanism.system
    set_data!(mechanism, data)
    set_solution!(mechanism, sol)
    setentries!(mechanism)
    return full_vector(system)
end

# function finitediff_data_jacobian(mechanism::Mechanism, data::AbstractVector,
#         sol::AbstractVector; δ = 1e-5, verbose = false)
#     FiniteDiff.finite_difference_jacobian(data -> evaluate_residual!(deepcopy(mechanism), data, sol), data)
# end

function finitediff_data_jacobian(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-5, verbose = false)
    mechanism = deepcopy(mechanism)
    Nd = data_dim(mechanism, attjac=false)
    N = length(mechanism)
    jac = zeros(N, Nd)

    setdata!(mechanism, deepcopy(data))
    setsolution!(mechanism, deepcopy(sol))

    for i = 1:Nd
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
