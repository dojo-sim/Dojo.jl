# Defining a Controller

Here, we explain how to write a controller and simulate its effect on a dynamical system, i.e., a [`Mechanism`](@ref).
We focus on the stabilization of the cartpole, which has two joints but only a single input on the cart. The controller is a method that always takes 2 input arguments:
- a [`Mechanism`](@ref),
- an integer `k` indicating the current simulation step.
For the cartpole, the controller computes the control input based on the current state `x`, the goal state `x_goal` and an LQR controller. The simulation step is not used in this example.

There are three ways to apply inputs to the system
- set an input directly to a joint
- set a set of inputs to all joints of the mechanism
- set an external force on bodies

```julia
# ### Setup
using Dojo
using DojoEnvironments

# ### Mechanism
mechanism = get_mechanism(:cartpole)

# ### Controller
K = [-0.948838; -2.54837; 48.6627; 10.871]

function controller!(mechanism, k)
    ## Target state
    x_goal = [0;0; 0.0;0]

    ## Current state
    x = get_minimal_state(mechanism)

    ## Control inputs
    u = -K' * (x - x_goal)

    ## 3 ways to set input:

    ## 1: get joint and set input
    cart_joint = get_joint(mechanism, :cart_joint)
    set_input!(cart_joint, [u])

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
```
