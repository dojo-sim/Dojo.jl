function step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::Vector{T}, u::Vector{T}; 
    opts=InteriorPointOptions{T}()) where {T,Nn,Ne,Nb,Ni}

    # set state
    setState!(mechanism, z)

    # set control
    setControl!(mechanism, u)

    # solve the 1-step simulation problem
    mehrotra!(mechanism, opts=opts)

    # extract the next state
    z̄ = getNextState(mechanism)
    return z̄
end