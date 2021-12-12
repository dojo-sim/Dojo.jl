################################################################################
# Hopper
################################################################################
struct Hopper end 

function hopper(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    control_scaling=Diagonal(ones(3)),
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer(),
    info=nothing,
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = gethopper(Δt=dt, g=g)
    initializehopper!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = 3
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e-3 * ones(nu)), high=(1.0e-3 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx) 
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = hopper_control_mask() 
    build_robot(vis, mechanism)

    TYPES = [Hopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
        
    return env
end

function hopper_control_mask()
    [0 0 0 1 0 0 0;
	 0 0 0 0 1 0 0;
	 0 0 0 0 0 0 1]
end

function hopper_nominal_max(; foot_radius=0.05)
    # initial state
    x1b1 = [0.0; 0.0; 0.5 + foot_radius]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; 0.0 + foot_radius]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function hopper_offset_max(x_shift, y_shift, z_shift)
    z = hopper_nominal_max()
    shift = [x_shift; y_shift; z_shift]
    z[1:3] += shift
    z[13 .+ (1:3)] += shift
    return z
end