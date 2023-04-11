# ### Setup
# PKG_SETUP
using Dojo
using DojoEnvironments

# ### List of all convenience mechanisms
mechanisms = [
    :ant, 
    :atlas,
    :block, 
    :block2d,
    :cartpole,
    :dzhanibekov,
    :fourbar, 
    :halfcheetah,
    :hopper, 
    :humanoid,
    :npendulum,
    :nslider,
    :panda,
    :pendulum,
    :quadruped,
    :raiberthopper,
    :slider,
    :snake,
    :sphere,
    :tippetop,
    :twister, 
    :walker,
    :youbot,
]

# ### Select mechanism
name = :ant

# ### Get mechanism (check DojoEnvironment/mechanisms files for kwargs)
mech = get_mechanism(name) 

# ### Initialize mechanism (check DojoEnvironment/mechanisms files for kwargs)
initialize!(mech, name)

# ### Simulate mechanism
storage = simulate!(mech, 5, record=true)
    
# ### Visualize mechanism
vis = visualize(mech, storage)
render(vis)
