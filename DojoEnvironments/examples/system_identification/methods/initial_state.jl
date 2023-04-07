function initial_state(model::Symbol; kwargs...)
    eval(Symbol(:initial_state_, model))(; kwargs...)
end

function initial_state_sphere(;
	xlims=[[0,0,0], [1,1,0.2]],
	vlims=[-1ones(3), 1ones(3)],
	ωlims=[-5ones(3), 5ones(3)],)
	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	return Dict(
		:position => x,
		:velocity => v ,
		:angular_velocity => ω)
end

function initial_state_block2d(;
		xlims=[[0,0.2], [1,0.4]],
        vlims=[-ones(2), ones(2)],
		θlims=[-π, π],
        ωlims=[-10, 10],)
	x = xlims[1] + rand(2) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(2) .* (vlims[2] - vlims[1])
	θ = θlims[1] + rand() * (θlims[2] - θlims[1])
	ω = ωlims[1] + rand() * (ωlims[2] - ωlims[1])
	return Dict(
		:position => x,
		:velocity => v ,
		:orientation => θ,
		:angular_velocity => ω)
end

function initial_state_block(;
		xlims=[[0,0,0.2], [1,1,0.4]],
        vlims=[-2ones(3), [2,2,-1.]],
        ωlims=[-6ones(3), 6ones(3)],)
	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	return Dict(
		:position => x,
		:orientation => rand(QuatRotation).q,
		:velocity => v ,
		:angular_velocity => ω)
end

function initial_state_nerf(;
	   nerf::Symbol=:bunny,
	   xlims=[[0,0,0.2], [1,1,0.4]],
       vlims=[[-2,-2,-0.5], [2,2,0.]],
       ωlims=[-6ones(3), 6ones(3)],)

	x = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	v = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	ω = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])

	return Dict(
		:position => x,
		:velocity => v ,
		:orientation => rand(QuatRotation).q,
		:angular_velocity => ω
		)
end

function initial_state_nerf_sphere(;
	r_sphere=1.0,
	v_sphere=1.0,
	nerf_position=[0,0.2,0.35],
	)

	q_nerf = Quaternion(normalize(rand(4))...)

	return Dict(
		:nerf_position => nerf_position,
		:sphere_position => r_sphere * [1,0,0.],
		:sphere_velocity => v_sphere * [-1,0,0.],
		:nerf_orientation => q_nerf,
		:sphere_orientation => one(Quaternion),
		)
end
