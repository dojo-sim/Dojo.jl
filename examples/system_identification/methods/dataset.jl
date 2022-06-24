################################################################################
# synthetic dataset
################################################################################

# Generate & Save Dataset
function generate_dataset(model::Symbol;
		N::Int=10,
		H=2.0,
		ctrl!::Function=(m,k) -> nothing,
		opts=SolverOptions(btol=1e-6, rtol=1e-6),
		init_kwargs=Dict(), # xlims, vlims, θlims, ωlims...
		mech_kwargs=Dict(), # timestep, gravity, friction_coefficient, radius, side...
		vis::Visualizer=Visualizer(),
		sleep_ratio = 0.0,
		show_contact = true,
		)
    mechanism = get_mechanism(model; mech_kwargs...);
    trajs = []
    for i = 1:N
		state = initial_state(model; init_kwargs...)
        initialize!(mechanism, model; state...)
        storage = simulate!(mechanism, H, ctrl!, record=true, opts=opts)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end

	data = get_data(mechanism)
    params = Dict(:N => N, :H => H, :timestep => timestep, :gravity => gravity, :data => data)
    jldsave(joinpath(@__DIR__, "..", "data", "synthetic", "dataset",
		datafilename(model; N=N, mech_kwargs...));
        params=params, trajs=trajs)
    return nothing
end

# Load Dataset
function open_dataset(model::Symbol; experiment_type="synthetic", kwargs...)
    dataset = jldopen(joinpath(@__DIR__, "..", "data", experiment_type, "dataset",
		datafilename(model; kwargs...)))
    params = dataset["params"]
    trajs = dataset["trajs"]
    JLD2.close(dataset)
    return params, trajs
end


################################################################################
# hardware dataset
################################################################################
function convert_csv_to()
	for k = 0:569
		csv_path = joinpath(module_dir(), "examples", "system_identification",
			"data", "hardware", "tosses_csv", "$(k).csv")
		jld2_path = joinpath(module_dir(), "examples", "system_identification",
			"data", "hardware", "tosses_jld2", "$(k).jld2")
		toss = split.(readlines(csv_path), ",")
		toss = [parse.(Float64, t) for t in toss]
		jldsave(jld2_path; toss=toss)
	end
end

# Convert the JLD2 files into trajectories and store them into a dataset
function raw_data_to_trajectory(toss, timestep::T; S=1) where T
	N = length(toss)
	M = Int((N-1 - (N-1)%S) / S)

	z = Vector{Vector{T}}()
	for i = 1:M
		vec1 = toss[1+S*(i-1)]
		vec2 = toss[1+S*i]

		x1 = vec1[1:3]
		q1 = vec1[4:7]

		x2 = vec2[1:3]
		q2 = vec2[4:7]

		v15 = (x2 - x1) / timestep
		ω15 = angular_velocity(Quaternion(q1...), Quaternion(q2...), timestep)

		z2 = [x2; v15; q2; ω15]
		push!(z, z2)
	end
	return z
end

function generate_hardware_dataset(model::Symbol;
		N::Int=10,
		H=1.0,
		opts=SolverOptions(btol=1e-6, rtol=1e-6),
		mech_kwargs=Dict(), # timestep, gravity, friction_coefficient, radius, side...
		vis::Visualizer=Visualizer(),
		sleep_ratio = 0.0,
		show_contact = true,
		S=1,
		)
	mechanism = get_mechanism(model; mech_kwargs...);
    trajs = []
    for i = 1:N
		file = jldopen(joinpath(module_dir(), "examples", "system_identification", "data",
			"hardware", "tosses_jld2", "$(i).jld2"))
		toss = file["toss"]
		z = raw_data_to_trajectory(toss, timestep, S=S)
		storage = generate_storage(mechanism, z)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end
	data = [
		0.2, 0, +1, +1, -1,
		0.2, 0, +1, -1, -1,
		0.2, 0, -1, +1, -1,
		0.2, 0, -1, -1, -1,
		0.2, 0, +1, +1, +1,
		0.2, 0, +1, -1, +1,
		0.2, 0, -1, +1, +1,
		0.2, 0, -1, -1, +1]
    params = Dict(:N => N, :H => H, :timestep => mechanism.timestep, :gravity => gravity, :data => data)
	jldsave(joinpath(@__DIR__, "..", "data", "hardware", "dataset",
		datafilename(model; N=N, mech_kwargs...));
		params=params, trajs=trajs)
    return nothing
end
