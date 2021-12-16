################################################################################
# Walker2d
################################################################################
struct Walker2d end

function walker2d(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    cf::T=1.9, spring=0.0, damper=0.1,
    s::Int=1, contact::Bool=true, info=nothing, vis::Visualizer=Visualizer(),
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = getwalker2d(Δt=dt, g=g, cf=cf, spring=spring, damper=damper, contact=contact)
    initializewalker2d!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max
        nx = maxCoordDim(mechanism)
    end
    nu = 6
    no = nx

    # values taken from Mujoco's model, combining the control range -1, 1 and the motor gears.
    aspace = BoxSpace(nu, low=(-ones(nu)), high=(ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(nu, 3) I(nu)]
    motor_gear = [100, 100, 100, 100, 100, 100.]
    control_scaling = Diagonal(dt * motor_gear)

    build_robot(vis, mechanism)

    TYPES = [Walker2d, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{Walker2d}; x=nothing, reset_noise_scale = 0.005)
    if x != nothing
        env.x .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :walker2d, z = 0.25)
        x0 = getMinState(env.mechanism)
        nx = minCoordDim(env.mechanism)

        low = -reset_noise_scale
        high = reset_noise_scale
        x = x0 + (high - low) .* rand(env.rng[1], nx) .+ low # we ignored the normal distribution on the velocities
        z = min2max(env.mechanism, x)
        setState!(env.mechanism, z)
        if env.mode == :min
            env.x .= getMinState(env.mechanism)
        elseif env.mode == :max
            env.x .= getMaxState(env.mechanism)
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end

function _get_obs(env::Environment{Walker2d}; full_state::Bool=false)
    full_state && (return env.x)
    nx = minCoordDim(env.mechanism)
    if env.mode == :min
        o = env.x
    elseif env.mode == :max
        o = max2min(env.mechanism, env.x)
    end
    # clamp velocities and remove x position
    ind = velocity_index(env.mechanism)
    vel = clamp.(o[ind], -10, 10)
    o[ind] .= vel
    o = o[[1;3:nx]]
    return o
end

function cost(env::Environment{Walker2d}, x, u;
        alive_bonus=0.1)

    if env.mode == :min
        x_velocity = -x[5]
    else
        i_torso = findfirst(body -> body.name == "torso", collect(env.mechanism.bodies))
        z_torso = x[(i_torso-1)*13 .+ (1:13)]
        x_velocity = -z_torso[4]
    end
    c = -x_velocity/10 # -forward velocity
    c -= alive_bonus
    c += 1e-4 * u'*u
    return c
end

function is_done(::Environment{Walker2d}, x)
    nx = minCoordDim(env.mechanism)
    if env.mode == :min
        x0 = x
    elseif env.mode == :max
        x0 = max2min(env.mechanism, x)
    end
    height = x0[1]
    ang = x0[3]
    done = !((height > 0.8) && (height < 2.0) && (abs(ang) < 1.0))
    return done
end


#
# env.x
# i_torso = findfirst(body -> body.name == "torso", collect(env.mechanism.bodies))
# z_torso = z[(i_torso-1)*13 .+ (1:13)]
# x_velocity = z_torso[4]
# z[3*13 + 4] = 324.0
# z
# setState!(env.mechanism, z)
#
# initialize!(env.mechanism, :walker2d, x = 111.0, z = 1.0, θ=0.18)
# x = getMinState(env.mechanism)
# z = getMaxState(env.mechanism)
# is_done(env, x)
#
