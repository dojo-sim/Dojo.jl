# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
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
	# u = 0.5 * mechanism.Δt * [szeros(6); sones(nu-6)]
	u = 0.3*[zeros(6); mech.Δt; zeros(nu-7)]
	setControl!(mech, u)
	return
end
Random.seed!(0)
mech = getmechanism(:humanoid, Δt=0.01, g=0.0, spring=0.0, damper=0.0, contact=false)
initialize!(mech, :humanoid)
ϵ = 1e-14
storage = simulate!(mech, 1.0, ctrl!, record=true, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))
visualize(mech, storage, vis=vis)

function astronaut_simulation(mech::Mechanism; tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.05)

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
		@show Δt[i]
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
		@show Δt[i]
		mech, storage, tcompute = astronaut_simulation(Nsim=Nsim, Δt=Δt[i],
			tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ, control_amplitude=control_amplitude)
		speed[i], energy[i] = process_energy_runs(mech, storage, tcompute, tctrl)
	end
	return speed, energy
end

################################################################################
# momentum
################################################################################

Δt = [0.10,]
mech, storage, tcompute = astronaut_simulation(;Nsim=1, Δt=1e-1, tsim=10.0, tctrl=1.0, seed=0, control_amplitude=0.05)
visualize(mech, storage[1], vis=vis)

Δt = [0.10, 0.03, 0.01, 0.003, 0.001]
speed_dj, ener_dj = benchmark_energy(Δt; Nsim=1, tsim=2.0, tctrl=1.0, seed=0, ϵ=1e-14)
speed_dj, plin_dj, pang_dj = benchmark_momentum(Δt; Nsim=1, tsim=1.0, seed=0, ϵ=1e-14)


# Saving results
path = joinpath(@__DIR__, "results/astronaut_momentum.jld2")
jldsave(path, speed=speed_dj, plin=plin_dj, pang=pang_dj)

# Recover results
file = jldopen(joinpath(@__DIR__, "../MuJoCoBenchmark.jl/results/astronaut_momentum.jld2"))
speed_mj, plin_mj, pang_mj, ener_mj = file["speed"], file["plin"], file["pang"], file["ener"]
file = jldopen(joinpath(@__DIR__, "results/astronaut_momentum.jld2"))
speed_dj, plin_dj, pang_dj = file["speed"], file["plin"], file["pang"]

plt = plot(layout=(1,2), legend=false, xlims=(1e-2,2e4), ylims=(1e-15,1e2))
plot!(plt[1,1], speed_dj, plin_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed_dj, pang_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed_dj, plin_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed_dj, pang_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)

plot!(plt[1,1], speed_mj, plin_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed_mj, pang_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed_mj, plin_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed_mj, pang_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)

# PGFPlot
using PGFPlotsX
using LaTeXStrings

mj_plin_plot = @pgf PlotInc({thick, "color=orange", "mark options={orange}"}, Table(speed_mj, plin_mj))
dj_plin_plot = @pgf PlotInc({thick, "color=cyan", "mark options={cyan}"}, Table(speed_dj, plin_dj))
mj_pang_plot = @pgf PlotInc({thick, "color=orange", "mark options={orange}"}, Table(speed_mj, pang_mj))
dj_pang_plot = @pgf PlotInc({thick, "color=cyan", "mark options={cyan}"}, Table(speed_dj, pang_dj))
mj_legend = LegendEntry(raw"MuJoCo")
dj_legend = LegendEntry(raw"Dojo")

axs1 = @pgf LogLogAxis({
		ymin=1e-15,
		ymax=1e1,
		title="linear momentum",
		xlabel=raw"$\times$ faster than real time",
		ylabel="N.s",
		},
		dj_plin_plot, dj_legend, mj_plin_plot, mj_legend,
		)

axs2 = @pgf LogLogAxis({
		ymin=1e-15,
		ymax=1e1,
		title="angular momentum",
		xshift=-20.0,
		# ticks="none",
		yticklabels={},
		xlabel=raw"$\times$ faster than real time",
		# ylabel="(N.s)",
		},
		dj_pang_plot, dj_legend, mj_pang_plot, mj_legend,
		)

gp = @pgf PGFPlotsX.GroupPlot(
    {group_style = {group_size="2 by 1"},
	legend_style = {nodes="{scale=0.70, transform shape}", anchor="east", at="{(0.95,0.35)}"},
    y_label_style={at="{(axis description cs:-0.00,.5)}", anchor="north"},
	# height = "6cm",
	# width = "6cm",
	# no_markers,
	# xlabel=raw"$\times$ faster than real time",
	# legend_pos="north west",
    },
    axs1, axs2)

filename = joinpath(@__DIR__, "figures/astronaut_momentum.tikz")
pgfsave(filename, gp; include_preamble = true, dpi = 150)




################################################################################
# energy
################################################################################

Δt = [0.10, 0.02, 0.005]
ener_traj = []
for i = 1:3
	mech, storage, tcompute = astronaut_simulation(;Nsim=1, Δt=Δt[i], tsim=100.0, tctrl=1.0, seed=0, control_amplitude=0.02)
	stride = Int(floor(1/Δt[i]))
	ener_t = [kineticEnergy(mech, storage[1], i) for i in stride+2:stride:length(storage[1])]
	push!(ener_traj, ener_t .- ener_t[1])
	@show i
end
ener_traj

# Saving results
path = joinpath(@__DIR__, "results/astronaut_energy.jld2")
jldsave(path, ener_traj=ener_traj)


# Recover results
file = jldopen(joinpath(@__DIR__, "../MuJoCoBenchmark.jl/results/astronaut_energy.jld2"))
ener_traj_mj = file["ener_traj"]
file = jldopen(joinpath(@__DIR__, "results/astronaut_energy.jld2"))
ener_traj_dj = file["ener_traj"]

plt = plot(layout=(3,1), legend=false, xlims=(0,100))
for i = 1:3
	plot!(plt[i,1],
		range(0,stop=100,length=length(ener_traj_mj[i])),
		ener_traj_mj[i], linewidth=3, markersize=6)
	plot!(plt[i,1],
		range(0,stop=100,length=length(ener_traj_dj[i])),
		ener_traj_dj[i], linewidth=3, markersize=6)
end
display(plt)



# PGFPlot
using PGFPlotsX
using LaTeXStrings

mj_ener_plot1 = @pgf PlotInc({very_thick, "color=orange", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_mj[1])),
		ener_traj_mj[1]))
mj_ener_plot2 = @pgf PlotInc({very_thick, "color=orange", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_mj[2])),
		ener_traj_mj[2]))
mj_ener_plot3 = @pgf PlotInc({very_thick, "color=orange", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_mj[3])),
		ener_traj_mj[3]))

dj_ener_plot1 = @pgf PlotInc({very_thick, "color=cyan", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_dj[1])),
		ener_traj_dj[1]))
dj_ener_plot2 = @pgf PlotInc({very_thick, "color=cyan", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_dj[2])),
		ener_traj_dj[2]))
dj_ener_plot3 = @pgf PlotInc({very_thick, "color=cyan", "mark options={orange}"},
	Table(range(0,stop=100,length=length(ener_traj_dj[3])),
		ener_traj_dj[3]))

mj_legend = LegendEntry(raw"MuJoCo")
dj_legend = LegendEntry(raw"Dojo")


axs1 = @pgf Axis({
		title=raw"simulation rate $10$ Hz",
		ylabel="energy J",
		xticklabels={},
		},
		dj_ener_plot1, dj_legend, mj_ener_plot1, mj_legend,
		)

axs2 = @pgf Axis({
		title=raw"simulation rate $100$ Hz",
		ylabel="energy J",
		yshift=13.0,
		xticklabels={},
		},
		dj_ener_plot2, dj_legend, mj_ener_plot2, mj_legend,
		)

axs3 = @pgf Axis({
		title=raw"simulation rate $1000$ Hz",
		xlabel="time s",
		ylabel="energy J",
		yshift=13.0,
		},
		dj_ener_plot3, dj_legend, mj_ener_plot3, mj_legend,
		)

gp = @pgf PGFPlotsX.GroupPlot(
    {group_style = {group_size="1 by 3"},
	legend_style = {nodes="{scale=0.70, transform shape}"},
	title_style={at="{(axis description cs:+0.5,0.95)}", anchor="south"},
	y_label_style={at="{(axis description cs:+0.05,0.5)}", anchor="north"},
    x_label_style={at="{(axis description cs:+0.5,0.15)}", anchor="north"},
	xmin=0,
	xmax=100,
	height = "3.5cm",
	width = "7.5cm",
	no_markers,
	# xlabel=raw"$\times$ faster than real time",
	legend_pos="north west",
    },
    axs1, axs2, axs3)

filename = joinpath(@__DIR__, "figures/astronaut_energy.tikz")
pgfsave(filename, gp; include_preamble = true, dpi = 150)
