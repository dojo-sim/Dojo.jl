# PREAMBLE

# PKG_SETUP

# ## Setup
using Dojo

# ## Mechanism
mechanism = get_mechanism(:pendulum, 
    timestep=0.01, 
    gravity=-9.81, 
    spring)

# ## Controller
function controller!(mechanism, k)
    ## Target state
    x_goal = [1.0 * π; 0.0]

    ## Current state
    x = get_minimal_state(mechanism) 

    ## Gains 
    K = [10.0 0.5] * 0.1

    off = 0
    for (i, eqc) in enumerate(mechanism.joints)
        nu = input_dimension(eqc)
        ## Get joint configuration + velocity
        xi = x[off .+ (1:2nu)]
        xi_goal = x_goal[off .+ (1:2nu)]
        
        ## Control
        ui = -K * (xi - xi_goal) 
        set_input!(eqc, ui)

        off += nu
    end
end

# ##Simulate
initialize!(mechanism, :pendulum, 
    ϕ1=0.0 * π, 
    ω1=0.0)

storage = simulate!(mechanism, 10.0, controller!, 
    record=true, 
    verbose=true)

# ## Visualize
vis = Visualizer()
render(vis)
visualize(mechanism, storage, 
    vis=vis)
