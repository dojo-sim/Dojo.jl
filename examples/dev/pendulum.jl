using Dojo 

# Open visualizer
vis=visualizer()
open(vis)
render(vis)

# Mechanism
mechanism = get_mechanism(:pendulum, timestep=0.01, gravity=0.0 * -9.81, spring=0.0)

# Controller
function controller!(mechanism, k)
    # target state
    x_goal = [1.0 * π; 0.0]

    # current state
    x = get_minimal_state(mechanism) 

    # gains 
    K = [10.0 0.5] * 0.1

    off = 0
    for (i, eqc) in enumerate(mechanism.joints)
        nu = input_dimension(eqc)
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
initialize!(mechanism, :pendulum, ϕ1 = 0.0 * π, ω1 = 0.0)

# Open visualizer
storage = simulate!(mechanism, 10.0, controller!, record=true, verbose=true)

# Visualize
visualize(mechanism, storage, vis=vis)
