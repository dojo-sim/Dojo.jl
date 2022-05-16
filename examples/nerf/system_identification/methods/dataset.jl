################################################################################
# Generate & Save Dataset
################################################################################

function generate_dataset(model::Symbol;
		N::Int=10, H=2.0, timestep=0.01, gravity=-9.81,
		opts=SolverOptions(btol=1e-6, rtol=1e-6),
		init_kwargs=Dict(), # xlims, vlims, θlims, ωlims...
		mech_kwargs=Dict(), # friction_coefficient, radius, side...
		vis::Visualizer=Visualizer(),
		sleep_ratio = 0.0,
		show_contact = true,
		)
    mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; mech_kwargs...);
    trajs = []
    for i = 1:N
		state = initial_state(model; init_kwargs...)
        initialize!(mechanism, model; state...)
        storage = simulate!(mechanism, H, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end

	data = get_data(mechanism)
    params = Dict(:N => N, :H => H, :timestep => timestep, :gravity => gravity, :data => data)
    jldsave(joinpath(@__DIR__, "..", "data", "dataset", datafilename(model; N = N, mech_kwargs...));
        params=params, trajs=trajs)
    return nothing
end

function generate_hardware_dataset(;N::Int=10,
		sleep_ratio = 0.0,
		show_contact = false,
		s::Int=1,
		scale=0.5,
		)
	H = 1.00
	timestep= 1/148 * s
	gravity_scaled = -9.81 * 20*scale

    mechanism = get_mechanism(:aprilcube, timestep=timestep, gravity=gravity_scaled, side=2scale);
    trajs = []
    for i = 1:N
		file = jldopen(joinpath(module_dir(), "examples", "system_identification", "data", "tosses_jld2", "$(i).jld2"))
		toss = file["toss"]
		z = toss_to_maximal_state(toss, timestep, s=s, scale=scale)
		storage = generate_storage(mech, z)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end
	data = get_data(mech) # TODO there is no ground-truth
    params = Dict(:N => N, :H => H, :timestep => timestep, :g => gravity_scaled, :data => data)
    jldsave(joinpath(@__DIR__, "..", "data", "dataset", datafilename(:aprilcube; N = N, s = s));
        params=params, trajs=trajs)
    return nothing
end

function toss_to_maximal_state(toss, timestep; s=1, scale=1.0)
	N = length(toss)
	M = Int((N-1 - (N-1)%s) / s)
	z = []
	for i = 1:M
		vec1 = toss[1+s*(i-1)]
		vec2 = toss[1+s*i]

		x1 = vec1[1:3] * scale
		q1 = vec1[4:7]

		x2 = vec2[1:3] * scale
		q2 = vec2[4:7]

		v15 = (x2 - x1) / timestep
		ϕ15 = angular_velocity(Quaternion(q1...),
			Quaternion(q2...), timestep)

		z2 = [x2; v15; q2; ϕ15]
		push!(z, z2)
	end
	return z
end

################################################################################
# Load Dataset
################################################################################
function open_dataset(model::Symbol; kwargs...)
    dataset = jldopen(joinpath(@__DIR__, "..", "data", "dataset", datafilename(model; kwargs...)))
    params = dataset["params"]
    trajs = dataset["trajs"]
    JLD2.close(dataset)
    return params, trajs
end
