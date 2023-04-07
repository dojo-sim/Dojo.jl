function process_momentum_runs(mechanism::Mechanism, storage::Vector{Storage{T,N}}, tcompute::Vector) where {T,N}
	H = mechanism.timestep * N
	Nsim = length(tcompute)
	speed = 0.0
	plin = 0.0
	pang = 0.0
	for i = 1:Nsim
		p = momentum(mechanism, storage[i], N)
		speed += H / tcompute[i]
		plin += norm(p[1:3])
		pang += norm(p[4:6])
	end
	speed /= Nsim
	plin /= Nsim
	pang /= Nsim
	return speed, plin, pang
end

function process_energy_runs(mechanism::Mechanism, storage::Vector{Storage{T,N}},
		tcompute::Vector, tctrl::T) where {T,N}
	H = mechanism.timestep * N
	Nsim = length(tcompute)
	speed = 0.0
	energy = 0.0
	for i = 1:Nsim
		speed += H / tcompute[i]
		initial = Int(floor(tctrl/mechanism.timestep+2))
		final = N
		energy += kinetic_energy(mechanism, storage[i], final) - kinetic_energy(mechanism, storage[i], initial)
	end
	speed /= Nsim
	energy /= Nsim
	return speed, energy
end