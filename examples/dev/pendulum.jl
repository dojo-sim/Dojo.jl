using Dojo 

# Open visualizer
vis = Visualizer()
open(vis)

# Mechanism
mechanism = get_mechanism(:pendulum, timestep=0.01, gravity=-9.81)

# Controller
function controller!(mechanism, k)
    # target state
    x_goal = [0.5 * π; 0.0]

    # current state
    x = get_minimal_state(mechanism) 

    # gains 
    K = [10.0 0.5]

    off = 0
    for (i, eqc) in enumerate(mechanism.joints)
        nu = control_dimension(eqc)
        # get joint configuration + velocity
        xi = x[off .+ (1:2nu)]
        xi_goal = x_goal[off .+ (1:2nu)]
        
        # control
        ui = -K * (xi - xi_goal) 
        set_input!(eqc, ui)

        off += nu
    end
end

@show [joint.name for joint in mechanism.joints]

# Simulate
initialize!(mechanism, :pendulum, ϕ1 = 0.0)
storage = simulate!(mechanism, 5.0, controller!, record=true, verbose=true)

# Visualize
visualize(mechanism, storage, vis=vis)

