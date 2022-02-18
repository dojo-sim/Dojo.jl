using Dojo 

# Open visualizer
vis = Visualizer()
open(vis)
render(vis)

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
        nu = control_dimension(eqc)
        # get joint configuration + velocity
        xi = x[off .+ (1:2nu)]
        xi_goal = x_goal[off .+ (1:2nu)]
        
        # control
        # ui = -K * (xi - xi_goal) 
        ui = 0.01 * sones(nu)
        set_input!(eqc, ui)

        off += nu
    end
end

# Mechanism
mechanism = get_mechanism(:pendulum, timestep=0.01, gravity=0.0 * -9.81, 
    spring=0.0,
    damper=0.0,
    qoffset=UnitQuaternion(RotZ(0.0 * π)))

# Simulate
initialize!(mechanism, :pendulum, ϕ1 = 0.5 * π, ω1 = 0.0)

# Open visualizer
storage = simulate!(mechanism, 10.0, controller!, record=true, verbose=true)

# Visualize
visualize(mechanism, storage, vis=vis)

z = get_maximal_state(storage)
x = [maximal_to_minimal(mechanism, zt) for zt in z]

plot(hcat(z...)[1:3,:]')
plot((hcat(x...)[1:2,:])')# yaxis=:log)