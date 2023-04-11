# ### Setup
# PKG_SETUP
using Dojo
using DojoEnvironments

# ### Mechanism
mechanism = get_mechanism(:pendulum)

# ### Controller
summed_error = 0 # for integral part

function controller!(mechanism, k)
    ## Target state
    x_goal = [π/2; 0.0]

    ## Current state
    x = get_minimal_state(mechanism)

    error = (x_goal-x)
    global summed_error += error[1]*mechanism.timestep

    ## Gains
    Kp = 25
    Ki = 50
    Kd = 5

    ## Control inputs
    u = Kp*error[1] + Ki*summed_error + Kd*error[2]
    set_input!(mechanism, [u])
end

# ### Simulate
initialize!(mechanism, :pendulum,
    angle=0.0, angular_velocity=0.0)

storage = simulate!(mechanism, 5.0, controller!, record=true)

# ### Visualize
vis = visualize(mechanism, storage)
render(vis)
