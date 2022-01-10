################################################################################
# Hopper
################################################################################
struct Cartpole end

function cartpole(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    s::Int=1, vis::Visualizer=Visualizer(), info=nothing, name::Symbol=:robot,
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = getcartpole(Δt=dt, g=g)
    initializecartpole!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max
        nx = maxCoordDim(mechanism)
    end
    nu = 1
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e-3 * ones(nu)), high=(1.0e-3 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [1.0 0.0]

    build_robot(vis, mechanism, name=name)

    TYPES = [Cartpole, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function cartpole_nominal_max( ;pendulum_length=1.0)
    # initial state
    x1b1 = [0.0; 0.0; 0.0]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; -0.5 * pendulum_length]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function cartpole_goal_max(;pendulum_length=1.0)
    # target state
    xTb1 = [0.0; 0.0; 0.0]
    vTb1 = [0.0; 0.0; 0.0]
    qTb1 = [1.0; 0.0; 0.0; 0.0]
    ωTb1 = [0.0; 0.0; 0.0]
    zTb1 = [xTb1; vTb1; qTb1; ωTb1]

    xTb2 = [0.0; 0.0; 0.5 * pendulum_length]
    vTb2 = [0.0; 0.0; 0.0]
    qTb2 = [0.0; 1.0; 0.0; 0.0]
    ωTb2 = [0.0; 0.0; 0.0]
    zTb2 = [xTb2; vTb2; qTb2; ωTb2]

    zT = [zTb1; zTb2]
end
