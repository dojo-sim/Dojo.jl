using OSFLoader
using Plots

results_dir = joinpath(OSFLoader.osf_loader_dir(), "assets/scan")

################################################################################
# bunny
################################################################################
bunny_nerf = OSFLoader.get_nerf_object(filename="bunny")

p = zeros(3)
xrange = p[1] .+ (-1:0.02:1)
yrange = p[2] .+ (-1:0.02:1)
zrange = p[3] .+ (-0.5:0.02:0.70)
z = 0.1
density = slice_density(bunny_nerf, xrange, yrange, z)

anim = nerf_scan(bunny_nerf, xrange, yrange, zrange, title="bunny density")
gif(anim, joinpath(results_dir, "bunny_scan.gif"), fps=15)


################################################################################
# bluesoap
################################################################################
bluesoap_nerf = OSFLoader.get_nerf_object(filename="bluesoap")

p = [-0.0207,  2.5801,  4.9728]
xrange = p[1] .+ (-3:0.05:3)
yrange = p[2] .+ (-3:0.05:3)
zrange = p[3] .+ (-3:0.05:3)
z = 0.1
density = slice_density(bluesoap_nerf, xrange, yrange, z)

anim = nerf_scan(bluesoap_nerf, xrange, yrange, zrange, title="bluesoap density")
gif(anim, joinpath(results_dir, "bluesoap_scan.gif"), fps=15)


################################################################################
# halfsoap
################################################################################
halfsoap_nerf = OSFLoader.get_nerf_object(filename="halfsoap")

p = [-2.5519,  2.4537,  4.8040]
xrange = p[1] .+ (-3:0.05:3)
yrange = p[2] .+ (-3:0.05:3)
zrange = p[3] .+ (-3:0.05:3)
z = 0.1
density = slice_density(halfsoap_nerf, xrange, yrange, z)

anim = nerf_scan(halfsoap_nerf, xrange, yrange, zrange, title="halfsoap density")
gif(anim, joinpath(results_dir, "halfsoap_scan.gif"), fps=15)
