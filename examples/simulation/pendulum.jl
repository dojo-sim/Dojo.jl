using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ## Setup
using Dojo

# ## Mechanism
mechanism = get_mechanism(:pendulum,
    timestep=0.01,
    gravity=-9.81,
    damper=5.0,
    spring=0.0)

# ## Controller
function controller!(mechanism, k)
    ## Target state
    x_goal = [1.0 * π; 0.0]

    ## Current state
    x = get_minimal_state(mechanism)

    ## Gains
    K = [5.0 0.5] * 0.1

    # Control inputs
    u =  -K * (x - x_goal)
    set_input!(mechanism, u)
end

# ##Simulate
initialize!(mechanism, :pendulum,
    angle=0.0 * π,
    angular_velocity=0.0);

storage = simulate!(mechanism, 2.0, controller!,
    record=true,
    verbose=true);

# ## Visualize
vis = Visualizer()
render(vis)
visualize(mechanism, storage, vis=vis);
