function astronaut_simulation(mech::Mechanism; 
		tsim=1.0, 
		tctrl=1.0, 
		seed=0, 
		ϵ=1.0e-14,
		control_amplitude=0.05)

    Random.seed!(seed)
    initialize!(mech, :humanoid)

	function ctrl!(mechanism, k)
		nu = input_dimension(mech)
		u = (k*mechanism.timestep < tctrl) * control_amplitude * [szeros(6); srand(nu-6)]
		set_input!(mech, u)
	    return
	end

    tcompute = @elapsed storage = simulate!(mech, tsim, ctrl!, 
		record=true,
		opts=SolverOptions(rtol=ϵ, btol=ϵ))

    return storage, tcompute
end

function astronaut_simulation(;
	Nsim=1, 
	timestep=1.0e-2, 
	gravity=0.0, 
	spring=0.0, 
	damper=0.0,
	tsim=1.0, 
	tctrl=1.0, 
	seed=0, 
	ϵ=1.0e-14, 
	control_amplitude=0.1)

    mech = get_mechanism(:humanoid, 
		timestep=timestep, 
		gravity=gravity, 
		springs=spring, 
		dampers=damper, 
		contact_feet=false)

	storage = []
	tcompute = zeros(Nsim)

	for i = 1:Nsim
		s, tcompute[i] = astronaut_simulation(mech; 
			tsim=tsim, 
			tctrl=tctrl, 
			seed=seed+i,
			ϵ=ϵ, 
			control_amplitude=control_amplitude)
		push!(storage, s)
	end
	return mech, [storage...], tcompute
end