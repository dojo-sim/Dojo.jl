# PREAMBLE

# PKG_SETUP

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

    off = 0
    for joint in mechanism.joints
        nu = input_dimension(joint)

        ## Get joint configuration + velocity
        xi = x[off .+ (1:2nu)]
        xi_goal = x_goal[off .+ (1:2nu)]

        ## Control
        ui = -K * (xi - xi_goal)
        set_input!(joint, ui)

        off += nu
    end
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
