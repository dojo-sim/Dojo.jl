################################################################################
# Hopper
################################################################################
struct HalfCheetah end 

function halfcheetah(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    cf::T=0.8, spring::T=0.0, damper::T=1.0, s::Int=1, contact::Bool=true, 
    info=nothing, vis::Visualizer=Visualizer(),
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = gethalfcheetah(Δt=dt, g=g, cf=cf, spring=spring, damper=damper)
    initializehalfcheetah!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = 6
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e-3 * ones(nu)), high=(1.0e-3 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx) 
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(6, 3) I(nu)]

    build_robot(vis, mechanism)

    TYPES = [HalfCheetah, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, control_mask,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
        
    return env
end
