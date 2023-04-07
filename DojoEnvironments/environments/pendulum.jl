"""
    Pendulum <: Environment

    classic system with one rotational degree of freedom
    https://underactuated.mit.edu/pend.html
"""
struct Pendulum end

function pendulum(; 
    representation=:minimal, 
    max_speed=8, 
    max_torque=8,
    timestep=0.05, 
    gravity=-10, 
    mass=1, 
    len=1, 
    springs=0,
    dampers=0, 
    seed=1, 
    vis=Visualizer(), 
    name=:robot,
    opts_step=SolverOptions(),
    opts_grad=SolverOptions(),
    T=Float64)

    mechanism = get_mechanism(:pendulum;
        timestep, 
        gravity, 
        mass, 
        len, 
        damper)

    initialize!(mechanism, :pendulum)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
        no = 3
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
        no = 13
    end
    nu = input_dimension(mechanism)

    high = [1, 1, max_speed]
    aspace = BoxSpace(input_dimension(mechanism), 
        low=[-timestep * max_torque], 
        high=[timestep * max_torque])
    ospace = BoxSpace(no, 
        low=-high, 
        high=high)
    rng = [MersenneTwister(seed),]

    x = Inf * ones(nx)
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = Inf * ones(nu)
    control_mask = ones(1, 1)
    control_scaling = Diagonal(ones(nu))
    build_robot(mechanism,
        vis=vis, 
        name=name)

    info = Dict(:maximal_speed => max_speed, :maximal_torque => max_torque)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{Pendulum, TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling, nx, nu, no, info,
        rng, vis, opts_step, opts_grad)
    return env
end

function Base.reset(env::Environment{Pendulum}; 
    x=nothing)

    initialize!(env.mechanism, :pendulum)

    if x != nothing
        env.state .= x
    else
        if env.representation == :minimal
            high = [π, 1]
            low = -high
            env.state .= rand(env.rng[1], env.num_states) .* (high .- low) .+ low
        elseif env.representation == :maximal
            env.state .= pendulum_nominal_max()
        end
        env.input_previous .= Inf
    end
    return get_observation(env)
end

function get_observation(env::Environment{Pendulum})
    if env.representation == :minimal
        θ, ω = env.state
        return [cos(θ), sin(θ), ω]
    else env.representation == :maximal
        return env.state
    end
end

function Base.step(env::Environment{Pendulum}, x, u; 
    gradients=false,
    attitude_decompress=false)
    mechanism = env.mechanism
    timestep= mechanism.timestep
    max_torque = env.info[:maximal_torque]

    x0 = x
    u0 = clamp.(u, -max_torque, max_torque)
    env.input_previous .= u0

    z0 = env.representation == :minimal ? minimal_to_maximal(mechanism, x0) : x0
    z1 = step!(mechanism, z0, timestep * u0; opts = env.opts_step)
    env.state .= env.representation == :minimal ? maximal_to_minimal(mechanism, z1) : z1

    # Compute cost function
    costs = cost(env, x0, u0)

    # Gradients
    if gradients
        if env.representation == :minimal
            fx, fu = get_minimal_gradients!(env.mechanism, z0, timestep * u0, 
                opts=env.opts_grad)
        elseif env.representation == :maximal
            fx, fu = get_maximal_gradients!(env.mechanism, z0, timestep * u0, 
                opts=env.opts_grad)
            if attitude_decompress
                A0 = attitude_jacobian(z0, length(env.mechanism.bodies))
                A1 = attitude_jacobian(z1, length(env.mechanism.bodies))
                fx = A1 * fx * A0'
                fu = A1 * fu
            end
        end
        env.dynamics_jacobian_state .= fx
        env.dynamics_jacobian_input .= timestep * fu
    end

    info = Dict()
    return get_observation(env), -costs, false, info
end

function angle_normalize(x)
    return ((x + 101π) % (2 * π)) - π
end

function pendulum_nominal_max()
    x1 = [0; 0; -0.5]
    v1 = [0; 0; 0]
    q1 = [1; 0; 0; 0]
    ω1 = [0; 0; 0]
    z1 = [x1; v1; q1; ω1]
end

function pendulum_goal_max()
    xT = [0; 0; 0.5]
    vT = [0; 0; 0]
    qT = [0; 1; 0; 0]
    ωT = [0; 0; 0]
    zT = [xT; vT; qT; ωT]
end

function cost(env::Environment{Pendulum}, x, u)
    if env.representation == :minimal
        θ, ω = x
        c = angle_normalize(θ - π)^2 + 1e-1 * ω^2 + 1e-3 * (u[1])^2 # angle_normalize enforces angle ∈ [-π, π]
        c = angle_normalize(θ - π)^2 + 1e-3 * ω^2 + 1e-3 * (env.mechanism.timestep * u[1])^2 # angle_normalize enforces angle ∈ [-π, π]
    else
        c = Inf
    end
    return c * 0.1
end

# function initialize_pendulum!(mechanism::Mechanism;
#     angle=0.7,
#     angular_velocity=0)
#     joint = mechanism.joints[1]
#     set_minimal_coordinates_velocities!(mechanism, joint;
#         xmin=[angle, angular_velocity])
# end