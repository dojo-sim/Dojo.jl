################################################################################
# build nerf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
nerf_object = py"generate_test_nerf"()

################################################################################
# Build Scan
################################################################################
xrange = -1.0:0.01:1.0
yrange = -1.0:0.01:1.0
zrange = -0.5:0.01:0.65
z = 0.1
density = slice_density(nerf_object, xrange, yrange, z)

plt = heatmap(xrange, yrange, density,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x", ylabel="y",
    title="bunny")

anim = @animate for z ∈ zrange
    points = slice_points(xrange, yrange, z)
    density = py"density_query"(nerf_object, points)
    density = reshape(density, length(xrange), length(yrange))
    plt = heatmap(xrange, yrange, density,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x", ylabel="y",
        title="bunny")
    display(plt)
end

nerf_dir = joinpath(module_dir(), "results")
gif(anim, joinpath(nerf_dir, "bunny_scan.gif"), fps=15)
