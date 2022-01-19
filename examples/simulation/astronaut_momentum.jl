using Dojo

# visualizer
vis = Visualizer()
open(vis)

################################################################################
# Astronaut
################################################################################
# humanoid
# multiple bodies
# no initial linear velocities
# no gravity
# no spring and damper
# random control
################################################################################

function ctrl!(mechanism, k)
	nu = controldim(mech)
	u = 0.5 * mechanism.Δt * [szeros(6); sones(nu-6)]
	setControl!(mech, u)
	return
end
Random.seed!(0)
mech = getmechanism(:humanoid, Δt=0.01, g=0.0, spring=0.0, damper=0.0, contact=false)
initialize!(mech, :humanoid)
ϵ = 1e-14
@profiler storage = simulate!(mech, 1.0, ctrl!, record=true, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))
visualize(mech, storage, vis=vis)

function astronaut_simulation(mech::Mechanism; tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.1)

    Random.seed!(seed)
    initialize!(mech, :humanoid)

	function ctrl!(mechanism, k)
		nu = controldim(mech)
		u = (k*mechanism.Δt < tctrl) * control_amplitude * mechanism.Δt * [szeros(6); srand(nu-6)]
		setControl!(mech, u)
	    return
	end
    tcompute = @elapsed storage = simulate!(mech, tsim, ctrl!, record=true,
		opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))
    return storage, tcompute
end

function astronaut_simulation(;Nsim::Int=1, Δt=1e-2, g=0.0, spring=0.0, damper=0.0,
		tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14, control_amplitude=0.1)
    mech = getmechanism(:humanoid, Δt=Δt, g=g, spring=spring, damper=damper, contact=false)
	storage = []
	tcompute = zeros(Nsim)
	for i = 1:Nsim
		s, tcompute[i] = astronaut_simulation(mech; tsim=tsim, tctrl=tctrl, seed=seed+i,
			ϵ=ϵ, control_amplitude=control_amplitude)
		push!(storage, s)
	end
	return mech, [storage...], tcompute
end

function process_momentum_runs(mechanism::Mechanism, storage::Vector{Storage{T,N}}, tcompute::Vector) where {T,N}
	H = mechanism.Δt * N
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
	H = mechanism.Δt * N
	Nsim = length(tcompute)
	speed = 0.0
	energy = 0.0
	for i = 1:Nsim
		speed += H / tcompute[i]
		initial = Int(floor(tctrl/mechanism.Δt+2))
		final = N
		energy += kineticEnergy(mechanism, storage[i], final) - kineticEnergy(mechanism, storage[i], initial)
	end
	speed /= Nsim
	energy /= Nsim
	return speed, energy
end

function benchmark_momentum(Δt::Vector; Nsim=1, tsim=1.0, tctrl=1.0, seed=0, ϵ=1e-14, control_amplitude=0.05)
	Nt = length(Δt)
	speed = zeros(Nt)
	plin = zeros(Nt)
	pang = zeros(Nt)
	for i = 1:Nt
		mech, storage, tcompute = astronaut_simulation(Nsim=Nsim, Δt=Δt[i],
			tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ, control_amplitude=control_amplitude)
		speed[i], plin[i], pang[i] = process_momentum_runs(mech, storage, tcompute)
	end
	return speed, plin, pang
end

function benchmark_energy(Δt::Vector; Nsim=1, tsim=2.0, tctrl=1.0, seed=0, ϵ=1e-14, control_amplitude=0.05)
	Nt = length(Δt)
	speed = zeros(Nt)
	energy = zeros(Nt)
	for i = 1:Nt
		mech, storage, tcompute = astronaut_simulation(Nsim=Nsim, Δt=Δt[i],
			tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ, control_amplitude=control_amplitude)
		speed[i], energy[i] = process_energy_runs(mech, storage, tcompute, tctrl)
	end
	return speed, energy
end

Δt = [0.10, 0.05, 0.01, 0.005, 0.001]
# Δt = [0.10,]
speed0, plin0, pang0 = benchmark_momentum(Δt; Nsim=1, tsim=1.0, seed=0, ϵ=1e-14)
speed0, ener0 = benchmark_energy(Δt; Nsim=128, tsim=2.0, tctrl=1.0, seed=0, ϵ=1e-14)
plt = plot(layout=(1,2), legend=false, xlims=(1e-2,2e4), ylims=(1e-15,1e2))
plot!(plt[1,1], speed0, plin0, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed0, pang0, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed0, plin0, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed0, pang0, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
