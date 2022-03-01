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
