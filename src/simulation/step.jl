function step!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::Vector{T}, u::Vector{T}; 
    opts=InteriorPointOptions{T}()) where {T,Nn,Ne,Nb,Ni}

    # set state
    set_state!(mechanism, z)

    # set control
    set_control!(mechanism, u)

    # solve the 1-step simulation problem
    mehrotra!(mechanism, opts=opts)

    # extract the next state
    z̄ = get_next_state(mechanism)
    return z̄
end