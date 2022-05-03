using Dojo
using OSFLoader
using PyCall

################################################################################
# bunny
################################################################################
bunny_nerf = OSFLoader.get_nerf_object(filename="bunny_trans")

p = zeros(3)
xrange = p[1] .+ (-1:0.02:1)
yrange = p[2] .+ (-1:0.02:1)
zrange = p[3] .+ (-0.5:0.02:0.70)
z = 0.1
density = slice_density(bluesoap_nerf, xrange, yrange, z)

heatmap(xrange, yrange, density,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x", ylabel="y",
    title="bunny")

anim = @animate for z ∈ zrange
    density = slice_density(bunny_nerf, xrange, yrange, z)
    plt = heatmap(xrange, yrange, density,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x", ylabel="y",
        title="bunny")
    display(plt)
end

results_dir = joinpath(Dojo.module_dir(), "examples/polytope/results")
gif(anim, joinpath(results_dir, "bunny_scan_.gif"), fps=15)


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

heatmap(xrange, yrange, density,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x", ylabel="y",
    title="bluesoap")

anim = @animate for z ∈ zrange
    density = slice_density(bluesoap_nerf, xrange, yrange, z)
    plt = heatmap(xrange, yrange, density,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x", ylabel="y",
        title="bluesoap")
    display(plt)
end

results_dir = joinpath(Dojo.module_dir(), "examples/polytope/results")
gif(anim, joinpath(results_dir, "bluesoap_scan_.gif"), fps=15)
