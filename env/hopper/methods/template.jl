function hopperState(; x::T=0.0, z::T=0.0, θ::T=0.0) where T
    mechanism = getmechanism(:hopper)
    initialize!(mechanism, :hopper, x=x, z=z, θ=θ)

    Nb = length(mechanism.bodies)
    x = zeros(13 * Nb)

    for (i, body) in enumerate(mechanism.bodies)
        x2 = body.state.x2[1]
        v15 = zeros(3)
        q2 = body.state.q2[1]
        ϕ15 = zeros(3)
        x[13 * (i-1) .+ (1:13)] = [x2;  v15; vector(q2); ϕ15]
    end
    return x
end
