################################################################################
# Convert the csv files gotten from ContactNets into the jld2 format
################################################################################

function convert_csv()
	for k = 0:569
		csv_path = joinpath(module_dir(), "examples", "real2sim", "data", "tosses_csv", "$(k).csv")
		jld2_path = joinpath(module_dir(), "examples", "real2sim", "data", "tosses_jld2", "$(k).jld2")
		toss = split.(readlines(csv_path), ",")
		toss = [parse.(Float64, t) for t in toss]
		jldsave(jld2_path;	toss=toss)
	end
end

# convert_csv()

################################################################################
# Convert the JLD2 files into trajectories and store them into a dataset
################################################################################

function toss2z(toss, Δt; S::Int=1)
	N = length(toss)
	M = Int((N-1 - (N-1)%S) / S)
	z = []
	for i = 1:M
		vec1 = toss[1+S*(i-1)]
		vec2 = toss[1+S*i]

		x1 = vec1[1:3]
		q1 = vec1[4:7]

		x2 = vec2[1:3]
		q2 = vec2[4:7]

		v15 = (x2 - x1) / Δt
		ϕ15 = ω_finite_difference(UnitQuaternion(q1...),
			UnitQuaternion(q2...), Δt)

		z2 = [x2; v15; q2; ϕ15]
		push!(z, z2)
	end
	return z
end

function build_pairs(z)
    pairs = []
	N = length(z)
    for t = 1:N-1
        z1 = z[t]
        z2 = z[t+1]
        pair = [z1, z2]
        push!(pairs, pair)
    end
    return pairs
end


function generate_hardware_dataset(;N::Int=10,
		sleep_ratio = 0.0,
		show_contact = true,
		S::Int=1,
		)
	H = 1.00
	Δt = 1/148 * S
	gscaled = -9.81*20

    mechanism = getmechanism(:box, Δt=Δt, g=gscaled);
    trajs = []
	pairs = []
    for i = 1:N
		file = jldopen(joinpath(module_dir(), "examples", "real2sim", "data", "tosses_jld2", "$(i).jld2"))
		toss = file["toss"]
		# z = toss2z(toss[1:7], Δt)
		# z = toss2z(toss[end-10:end], Δt)
		z = toss2z(toss, Δt, S=S)
		storage = generate_storage(mech, z)
		push!(pairs, build_pairs(z)...)
        push!(trajs, storage)
        visualize(mechanism, storage, vis=vis, show_contact=show_contact)
		sleep(H*sleep_ratio)
    end
	data = [
		0.2, 0,0,0, +1, +1, -1,
		0.2, 0,0,0, +1, -1, -1,
		0.2, 0,0,0, -1, +1, -1,
		0.2, 0,0,0, -1, -1, -1,
		0.2, 0,0,0, +1, +1, +1,
		0.2, 0,0,0, +1, -1, +1,
		0.2, 0,0,0, -1, +1, +1,
		0.2, 0,0,0, -1, -1, +1]
    params = Dict(:N => N, :H => H, :Δt => Δt, :g => gscaled, :data => data)
    jldsave(joinpath(@__DIR__, "dataset", filename(:hardwarebox; N = N, S = S));
        params=params, trajs=trajs, pairs=pairs)
    return nothing
end
