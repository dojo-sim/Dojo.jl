function step!(mechanism::Mechanism{T}, z::Vector{T}, u::Vector{T}; opts=SolverOptions{T}()) where T
    # set state
    set_state!(mechanism, z)

    # set control
    set_control!(mechanism, u)

    # solve the 1-step simulation problem
    mehrotra!(mechanism, opts=opts)

    # extract the next state
    z_next = get_next_state(mechanism)
    
    return z_next
end