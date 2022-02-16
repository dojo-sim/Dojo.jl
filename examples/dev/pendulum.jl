using Dojo 

# Open visualizer
vis = Visualizer()
open(vis)

# Mechanism
mechanism = get_mechanism(:pendulum, timestep=0.1, gravity=0.0 * -9.81)

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
initialize!(mechanism, :pendulum, ϕ1 = 0.0, ω1 = 1.0)
storage = simulate!(mechanism, 20.0, record=true, verbose=true)

z7 = [zt[7] for zt in z]
z = get_maximal_state(storage)
x = [maximal_to_minimal(mechanism, zt) for zt in z]

θ = Float64[]
for t = 1:length(z)
    push!(θ, nullspace_mask(mechanism.joints[1].rotational) * displacement(mechanism.joints[1].rotational, 
        mechanism.origin.state.x2[1], mechanism.origin.state.q2[1], 
        z[t][1:3], UnitQuaternion(z[t][6 .+ (1:4)]..., false), vmat=true))
end

plot(hcat(x...)')
plot!(hcat(θ...)')

# Visualize
visualize(mechanism, storage, vis=vis)

θ = []
for t = 1:length(z)
    push!(θ, vector(displacement(mechanism.joints[1].rotational, 
        mechanism.origin.state.x2[1], mechanism.origin.state.q2[1], 
        z[t][1:3], UnitQuaternion(z[t][6 .+ (1:4)]..., false), vmat=false)))
end

plot(hcat(θ...)')