# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## Setup
using Dojo
using DojoEnvironments
using ControlSystemsBase
using LinearAlgebra

# ## Mechanism
mechanism = get_mechanism(:cartpole)

# ## Controller
x0 = zeros(4)
u0 = zeros(2)
A, B = get_minimal_gradients!(mechanism, x0, u0)
Q = I(4)
R = I(1)
K = lqr(Discrete,A,B[:,1],Q,R)

function controller!(mechanism, k)
    ## Target state
    x_goal = [0;0; 0.0;0]

    ## Current state
    x = get_minimal_state(mechanism)

    ## Control inputs
    u = -K * (x - x_goal)
    set_input!(mechanism, [u;0]) # input only to cart
end

initialize!(mechanism, :cartpole; position=0, orientation=pi/4)
            
storage = simulate!(mechanism, 10.0, controller!, record=true)

    
# ## Visualize
visualize(mechanism, storage)
