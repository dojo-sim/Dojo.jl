function halfcheetahState(; x::T=0.0, z::T=0.0, θ::T=0.0) where T
    mechanism = get_mechanism(:halfcheetah)
    initialize!(mechanism, :halfcheetah, x=x, z=z, θ=θ)

    Nb = length(mechanism.bodies)
    x = zeros(13 * Nb)
    
    for (i, body) in enumerate(mechanism.bodies)
        x2 = body.state.x2
        v15 = zeros(3)
        q2 = body.state.q2
        ω15 = zeros(3)
        x[13 * (i-1) .+ (1:13)] = [x2;  v15; vector(q2); ω15]
    end
    return x
end