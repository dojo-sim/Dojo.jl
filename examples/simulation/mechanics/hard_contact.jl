# ### Setup
# PKG_SETUP_2
using Dojo
using DojoEnvironments
using Plots

# ### Simulate block
mech = get_mechanism(:block)
storage = simulate!(mech, 5, record=true)
vis = visualize(mech, storage)
render(vis)

# ### Contact interpenetration
distances = get_sdf(mech, storage) # distance from floor to each contact
minimum(minimum(distances)) # minimum distance of any corner to the ground

# ### Exemplary plot of one corner
plot(distances[1]) 