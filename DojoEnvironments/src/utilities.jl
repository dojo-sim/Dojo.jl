################################################################################
# Visuals
################################################################################
function visualize(env::Environment, traj::Vector{Vector{T}}; build::Bool=true) where T
	@assert size(traj[1]) == size(env.state)
    storage = generate_storage(env.mechanism, [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj])
    Dojo.visualize(env.mechanism, storage, vis=env.vis, build=build)
end

################################################################################
# Miscellaneous
################################################################################
type2symbol(H) = Symbol(lowercase(String(H.name.name)))

function get_control_mask(n_inputs, indices)
    n_outputs = length(indices)
    m = zeros(n_outputs, n_inputs)
    for (i,ind) in enumerate(indices)
        m[i,ind] = 1.0
    end
    return m
end
