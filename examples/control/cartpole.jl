# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## Setup
using Dojo
using DojoEnvironments
using ControlSystemsBase
using LinearAlgebra

# ## Mechanism
mechanism = get_mechanism(:cartpole,
    timestep=0.01,
    gravity=-9.81,
    damper=0.0,
    spring=0.0)

# ##Simulate
initialize!(mechanism, :cartpole,
    mode=:up);

set_maximal_configurations!(mechanism.bodies[1], mechanism.bodies[2], 
        Δx=[0.0; 0; 0.5], 
        Δq=Dojo.RotX(0))

# ## Controller
A, B = get_minimal_gradients!(mechanism)
Q = I(4)
R = I(1)
K = lqr(Discrete,A,B[:,1],Q,R)

function controller!(mechanism, t)
    ## Target state
    x_goal = [0;0; 0.0;0]

    ## Current state
    x = get_minimal_state(mechanism)

    # Control inputs
    u = -K * (x - x_goal)
    set_input!(mechanism, [u;0])
end

set_maximal_configurations!(mechanism.bodies[1], mechanism.bodies[2], 
            Δx=[0.0; 0; 0.5], 
            Δq=Dojo.RotX(0.1))
            
storage = simulate!(mechanism, 10.0, controller!,
    record=true,
    verbose=true);

    
# ## Visualize
vis = Visualizer()
render(vis)
visualize(mechanism, storage, 
    vis=vis);
