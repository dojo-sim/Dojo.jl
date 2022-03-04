# Defining an Environment

An [`Environment`](@ref) is a convienient object for applications like reinforcement learning and trajectory optimization. 

To demonstrate, we create the [`Dojo.Ant`](@ref) environment. First, we load (or [create](define_mechanism.md)) a mechanism:

```julia 
mechanism = get_mechanism(:ant) 
```

Next, we create an environment's other attributes.

Dimensions:
```julia
# set state dimension based on representation
if representation == :minimal
    nx = minimal_dimension(mechanism)
elseif representation == :maximal
    nx = maximal_dimension(mechanism)
end
# set control dimension
nu = 8
# set observation dimension
no = nx
```

Space (for limiting controls and observations):
```julia
# limit controls to [-1.0, 1.0]
aspace = BoxSpace(nu, 
    low=(-ones(nu)), 
    high=(ones(nu)))
# no limits on observations
ospace = BoxSpace(no, 
    low=(-Inf * ones(no)), 
    high=(Inf * ones(no)))
```

Random number:
```julia
rng = MersenneTwister(seed)
```

Dynamics data:
```julia
# state vector
z = get_maximal_state(mechanism)
x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z
# dynamics Jacobians
fx = zeros(nx, nx)
fu = zeros(nx, nu)
```

Control data: 
```julia
# control vector (previous)
u_prev = zeros(nu)
# control map transforms inputs from control to dynamics space
control_mask = [zeros(8, 6) I(nu)]
control_scaling = Diagonal(timestep * 150.0 * ones(nu))
control_map = control_mask' * control_scaling
```

Visuals: 
```julia 
# create a visualizer
vis = Visualizer() 
```

Solver options: 
```julia
# simulation options 
opts_step = SolverOptions()
# gradient options
opts_grad = SolverOptions() 
```

Environment:
```julia
TYPES = [Ant, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
env = Environment{TYPES...}(
    mechanism, 
    representation, 
    aspace, ospace,
    x, fx, fu,
    u_prev, 
    control_map,
    nx, nu, no,
    info,
    [rng], 
    vis,
    opts_sim, opts_grad)
```

With the environment instantiated, we can interact with it by overloading the following methods: 

Simulate environment forward one time step:
```julia
function step(env::Environment{Ant}, x, u; 
    gradients=false,
    attitude_decompress=false)

    # mechanism
    mechanism = env.mechanism

    # timestep 
    timestep = mechanism.timestep

    # copy current state
    x0 = copy(x)

    # cache current control
    env.input_previous .= u  # for rendering in Gym
	u_scaled = env.control_map * u

    # representation conversion
    z0 = env.representation == :minimal ? minimal_to_maximal(mechanism, x0) : x0

    # simulate one step
    z1 = step!(mechanism, z0, u_scaled; opts=env.opts_step)

    # representation conversion
    env.state .= env.representation == :minimal ? maximal_to_minimal(mechanism, z1) : z1

    # cost/reward
    reward = cost(env, z1, u_scaled)

    # check for done
    done = is_done(env, z1, u_scaled)

    # gradients
    if gradients
        if env.representation == :minimal
            fx, fu = get_minimal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        elseif env.representation == :maximal
            fx, fu = get_maximal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
            if attitude_decompress
                A0 = attitude_jacobian(z0, length(env.mechanism.bodies))
                A1 = attitude_jacobian(z1, length(env.mechanism.bodies))
                fx = A1 * fx * A0'
                fu = A1 * fu
            end
        end
        env.dynamics_jacobian_state .= fx
        env.dynamics_jacobian_input .= fu * env.control_map
    end

    # information
    info = Dict()

    return get_observation(env), reward, done, info
end
```

Return environment to nominal state:
```julia
function reset(env::Environment{Ant}; 
    x=nothing)

    # initialize
    initialize!(env.mechanism, type2symbol(Ant))

    if x != nothing
        env.state .= x
    else
        x = get_minimal_state(env.mechanism)
        if env.representation == :minimal
            set_maximal_state!(env.mechanism, minimal_to_maximal(env.mechanism, x))
            env.state .= x
        elseif env.representation == :maximal
            z = minimal_to_maximal(env.mechanism, x)
            set_maximal_state!(env.mechanism, z)
            env.state .= z
        end
        env.input_previous .= 0.0
    end

    return get_observation(env)
end
```

Observation for current environment state:
```julia
function get_observation(env::Environment{Ant})
    contact_force = Float64[]
    for contact in env.mechanism.contacts
        push!(contact_force, max(-1.0, min(1.0, contact.impulses[2][1])))
    end
    # include contact forces with state for observation
    return [env.state; contact_force]
end
```

Cost/reward associated with simulation step:
```julia 
function cost(env::Environment{Ant}, x, u)
    # forward reward
    v = x[4] # x-direction velocity
    forward_reward = 2.0 * v

    # control cost
	ctrl_cost = (0.05 * u' * u)[1]

    # contact cost
    contact_cost = 0.0

    for contact in mechanism.contacts
        contact_cost += 0.5 * 1.0e-3 * max(-1.0, min(1.0, contact.impulses[2][1]))^2.0
    end

	# survive_reward = 1.0
    survive_reward = 0.05

    # total reward
    reward = forward_reward - ctrl_cost - contact_cost + survive_reward
end
```

Determine if simulation should terminate:
```julia 
function is_done(env::Environment{Ant}, x) 
    !(all(isfinite.(env.state)) && (env.state[3] >= 0.2) && (env.state[3] <= 1.0))
end
```

### Random controls

We apply random controls to the robot via the environment interface:
```julia
y = [copy(env.state)] # state trajectory
for t = 1:100
    step(env, env.state, randn(env.num_inputs))
    push!(y, copy(env.state)) 
end
visualize(env, y)
```

The result should be something like this:
```@raw html
<img src="./assets/animations/ant_random.gif" width="300"/>
```
