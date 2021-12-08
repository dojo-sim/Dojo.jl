function step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
    ϵ::T = 1e-6, newtonIter::Int = 100, lineIter::Int = 10, verbose::Bool = true,
    btol::T = ϵ, undercut::T = Inf, ctrl!::Any = (m) -> nothing) where {T,Nn,Ne,Nb,Ni}

    # set the initial conditions: x1, v15, x2...
        # set x2, v15, q2, ϕ15
        # set x1, q1
        # set F2 and τ2 to zero
        # warm-start the solver
    setState!(mechanism, z)

    # set the controls in the equality constraints
        # apply the controls to each body's state
    setControl!(mechanism, u)

    # Apply a control policy, this is used to do puppet master control
    ctrl!(mechanism)

    # solve the 1 step simulation problem
    mehrotra!(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, verbose = verbose,
        opts=InteriorPointOptions(rtol=ϵ, max_iter=newtonIter, btol=btol, undercut=undercut, verbose=verbose))

    # extract the next state
    z̄ = getNextState(mechanism)
    return z̄
end