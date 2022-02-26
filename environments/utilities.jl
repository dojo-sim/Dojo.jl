################################################################################
# Visuals
################################################################################
function visualize(env::Environment, traj::Vector{Vector{T}}) where T
	@assert size(traj[1]) == size(env.state)
    storage = generate_storage(env.mechanism, [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj])
    visualize(env.mechanism, storage, 
        vis=env.vis)
end

################################################################################
# Miscellaneous
################################################################################
type2symbol(H) = Symbol(lowercase(String(H.name.name)))
