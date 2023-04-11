using Dojo
using DojoEnvironments
using JLD2

function convert_csv_to_jld2()
	timestep = 1/148
	mechanism = get_mechanism(:block; timestep)

	storages = []
	for k = 0:569
		csv_path = joinpath(@__DIR__, "tosses_csv", "$(k).csv")
		toss = split.(readlines(csv_path), ",")
		toss = [parse.(Float64, t) for t in toss]

		N = length(toss)
		z = Vector{Vector{Float64}}()
		for i=1:N-1
			z_i = toss[1+(i-1)]
			z_ip1 = toss[1+i]

			x1 = z_i[1:3]
			q1 = z_i[4:7]

			x2 = z_ip1[1:3]
			q2 = z_ip1[4:7]

			v15 = (x2 - x1) / timestep
			ω15 = Dojo.angular_velocity(Quaternion(q1...), Quaternion(q2...), timestep)

			push!(z, [x2; v15; q2; ω15])
		end

		push!(storages, generate_storage(mechanism, z))
	end

	jldsave(joinpath(@__DIR__, "..", "datasets", "real_block.jld2"); storages)

	return
end

convert_csv_to_jld2()