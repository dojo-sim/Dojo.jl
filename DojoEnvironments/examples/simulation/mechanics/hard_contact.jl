# ## Setup
using Dojo
using DojoEnvironments

# ## Simulate block
mech = get_mechanism(:block)
storage = simulate!(mech, 5, record=true)
visualize(mech, storage)

# ## Contact interpenetration
distances = get_sdf(mech, storage) # distance from floor to each contact
minimum(minimum(distances)) # minimum distance of any corner to the ground
plot(distances[1]) # exemplary plot of one corner