"""
    HalfCheetah <: Environment

    planar cheetah-like robot, based on https://gym.openai.com/envs/HalfCheetah-v2/
"""
struct HalfCheetah end

function halfcheetah(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=-9.81,
    friction_coefficient=0.4, 
    spring=[240, 180, 120, 180, 120, 60], 
    dampers=0,
    parse_dampers=true,
    limits=true,
    seed=1, 
    contact_feet=true, 
    contact_body=true,
    info=nothing, 
    vis=Visualizer(), 
    name=:robot,
    opts_step=SolverOptions(), 
    opts_grad=SolverOptions(),
    T=Float64)

    mechanism = get_halfcheetah(;
        timestep, 
        gravity, 
        friction_coefficient, 
        spring, 
        damper, 
        parse_damper, 
        contact_feet, 
        contact_body, 
        limits)
    initialize_halfcheetah!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = 6
    no = nx

    # values taken from Mujoco's model, combining the control range -1, 1 and the motor gears.
    aspace = BoxSpace(nu, 
        low=(-ones(nu)), 
        high=(ones(nu)))
    ospace = BoxSpace(no, 
        low=(-Inf * ones(no)), 
        high=(Inf * ones(no)))

    rng = MersenneTwister(seed)

    z = get_maximal_state(mechanism)
    x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(nu, 3) I(nu)]
    motor_gear = [120, 90, 60, 120, 60, 30.]
    control_scaling = Diagonal(timestep * motor_gear)

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [HalfCheetah, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function Base.reset(env::Environment{HalfCheetah}; 
    x=nothing, reset_noise_scale=0.1)
    if x != nothing
        env.state .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :halfcheetah, 
            body_position=[0, 0.25])
        x0 = get_minimal_state(env.mechanism)
        nx = minimal_dimension(env.mechanism)
        nz = maximal_dimension(env.mechanism)

        low = -reset_noise_scale
        high = reset_noise_scale
        x = x0 + (high - low) .* rand(env.rng[1], nx) .+ low # we ignored the normal distribution on the velocities
        z = minimal_to_maximal(env.mechanism, x)
        set_maximal_state!(env.mechanism, z)
        if env.representation == :minimal
            env.state .= get_minimal_state(env.mechanism)
        elseif env.representation == :maximal
            env.state .= get_maximal_state(env.mechanism)
        end
        env.input_previous .= 0.0
    end
    return get_observation(env)
end

function cost(env::Environment{HalfCheetah}, x, u;
        forward_reward_weight=1, 
        ctrl_cost_weight=0.1)

    if env.representation == :minimal
        x_velocity = -x[5]
    else
        i_torso = findfirst(body -> body.name == "torso", env.mechanism.bodies)
        z_torso = x[(i_torso-1) * 13 .+ (1:13)]
        x_velocity = -z_torso[4]
    end
    # @show mean(abs.(u))
    c = ctrl_cost_weight * u' * u - x_velocity * forward_reward_weight
    return c
end



# function initialize_halfcheetah!(mechanism::Mechanism{T}; 
#     body_position=[0, 0],  
#     body_orientation=0) where T

#     set_minimal_coordinates!(mechanism,
#                  get_joint(mechanism, :floating_joint),
#                  [body_position[2] + 0.576509, -body_position[1], -body_orientation + 0.02792])
#     for joint in mechanism.joints
#         (joint.name != :floating_joint) && set_minimal_coordinates!(mechanism, joint, zeros(input_dimension(joint)))
#     end
#     zero_velocity!(mechanism)
# end

# function halfcheetahState(; x::T=0, z::T=0, θ::T=0) where T
#     mechanism = get_mechanism(:halfcheetah)
#     initialize!(mechanism, :halfcheetah, x=x, z=z, θ=θ)

#     Nb = length(mechanism.bodies)
#     x = zeros(13 * Nb)
    
#     for (i, body) in enumerate(mechanism.bodies)
#         x2 = body.state.x2
#         v15 = zeros(3)
#         q2 = body.state.q2
#         ω15 = zeros(3)
#         x[13 * (i-1) .+ (1:13)] = [x2;  v15; vector(q2); ω15]
#     end
#     return x
# end


