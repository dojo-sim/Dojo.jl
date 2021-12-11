################################################################################
# Ant
################################################################################
struct Ant end 

function ant(; mode::Symbol=:min, dt::T=0.01, g::T=-9.81,
    cf::T=0.5, spring::T=0.0, damper::T=1.0, s::Int=1, 
    contact::Bool=true, contact_body=true,
    info=nothing, vis::Visualizer=Visualizer(),
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = getant(Δt=dt, g=g, cf=cf, spring=spring, damper=damper,contact=contact, contact_body=contact_body)
    initializeant!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = 8
    no = nx

    aspace = BoxSpace(nu, low=(-150 * dt * ones(nu)), high=(150 * dt * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx) 
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(8, 6) I(nu)]

    build_robot(vis, mechanism)

    TYPES = [Ant, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, control_mask,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
    return env
end
