function get_solution(mechanism::Mechanism{T}) where {T}
    sol = T[]
    for (i, joints) in enumerate(mechanism.joints)
        λ = joints.impulses[2]
        push!(sol, λ...)
    end
    for (i, body) in enumerate(mechanism.bodies)
        v25 = body.state.vsol[2]
        ϕ25 = body.state.ϕsol[2]
        push!(sol, [v25; ϕ25]...)
    end
    for (i, contact) in enumerate(mechanism.contacts)
        sol = get_solution(sol, contact)
    end
    return sol
end

function set_solution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
    off = 0
    for (i,joints) in enumerate(mechanism.joints)
        nλ = length(joints)
        λ = sol[off .+ (1:nλ)]; off += nλ
        joints.impulses[2] = λ
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
        sol, off = set_solution!(sol, off, contact)
    end
    return nothing
end

function get_solution(sol, contact::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    push!(sol, [s; γ]...)
end

function get_solution(sol, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    γ = contact.impulses[2]
    push!(sol, γ...)
end

function set_solution!(sol, off::Int, contact::RigidContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    s = sol[off .+ (1:N½)]; off += N½
    γ = sol[off .+ (1:N½)]; off += N½
    contact.impulses_dual[2] = s
    contact.impulses[2] = γ
    return sol, off
end

function set_solution!(sol, off::Int, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    γ = sol[off .+ (1:N)]; off += N
    contact.impulses[2] = γ
    return sol, off
end

function evaluate_residual!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
    system = mechanism.system
    set_data!(mechanism, data)
    set_solution!(mechanism, sol)
    set_entries!(mechanism)
    return full_vector(system)
end

function residual_dimension(mechanism::Mechanism)
    return sum(Vector{Int}(length.(mechanism.joints))) +
    	sum(Vector{Int}(length.(mechanism.bodies))) +
    	sum(Vector{Int}(length.(mechanism.contacts)))
end

function finite_difference_data_jacobian(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-5, verbose=false)
    mechanism = deepcopy(mechanism)
    for contact in mechanism.contacts
        initialize!(mechanism, contact)
    end
    Nd = data_dim(mechanism, attjac=false)
    Nr = residual_dimension(mechanism)
    jac = zeros(Nr, Nd)
    set_data!(mechanism, deepcopy(data))
    set_solution!(mechanism, deepcopy(sol))

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
