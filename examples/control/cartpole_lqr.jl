# ### Setup
# PKG_SETUP
using Dojo
using DojoEnvironments
using ControlSystemsBase
using LinearAlgebra

# ### Mechanism
mechanism = get_mechanism(:cartpole)

# ### Controller
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

    ## 3 ways to set input:

    ## 1: get joint and set input
    cart_joint = get_joint(mechanism, :cart_joint)
    set_input!(cart_joint, u)

    ## 2: set input for all joints at once
    ## set_input!(mechanism, [u;0]) # need to know joint order

    ## 3: direct external force on body
    ## cart = get_body(mechanism, :cart)
    ## set_external_force!(cart; force=[0;u;0])

end

# ### Simulate 
initialize!(mechanism, :cartpole; position=0, orientation=pi/4)
            
storage = simulate!(mechanism, 10.0, controller!, record=true)

    
# ### Visualize
vis = visualize(mechanism, storage)
render(vis)
