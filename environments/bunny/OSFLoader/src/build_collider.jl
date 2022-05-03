using Dojo
using OSFLoader
using PyCall

# run Dojo's activate.jl
# run using Dojo
# run OSFLoader's activate.jl
# run using OSFLoader, PyCall
# rerun Dojo's activate.jl
# this script needs to be executed in the Dojo module and Dojo Env after loading both Dojo and the OSFLOader modules

################################################################################
# build NeRf
################################################################################
bunny_nerf = OSFLoader.get_nerf_object(filename="bunny_trans")
halfsoap_nerf = OSFLoader.get_nerf_object(filename="halfsoap")
bluesoap_nerf = OSFLoader.get_nerf_object(filename="bluesoap")


# Nerf data
bunny_folder = joinpath(Dojo.module_dir(), "environments/bunny/deps")
halfsoap_folder = joinpath(Dojo.module_dir(), "environments/halfsoap/deps")
bluesoap_folder = joinpath(Dojo.module_dir(), "environments/bluesoap/deps")
outer_mesh = jldopen(joinpath(bunny_folder, "bunny_outer_mesh.jld2"))["mesh"]
inner_mesh = jldopen(joinpath(bunny_folder, "bunny_inner_mesh.jld2"))["mesh"]

################################################################################
# Offline
################################################################################

bunny = Dojo.SoftCollider(bunny_nerf, N=5000)
jldsave(joinpath(bunny_folder, "bunny_collider.jld2"), collider=bunny)

halfsoap = Dojo.SoftCollider(halfsoap_nerf, N=5000)
jldsave(joinpath(halfsoap_folder, "halfsoap_collider.jld2"), collider=halfsoap)

bluesoap = Dojo.SoftCollider(bluesoap_nerf, N=5000)
jldsave(joinpath(bluesoap_folder, "bluesoap_collider.jld2"), collider=bluesoap)


particles = zeros(Float32, 10, 3)
masses = OSFLoader.py"density_query"(halfsoap_nerf, particles)




################################################################################

using Plots

p = [-0.0207,  2.5801,  4.9728]
xrange = p[1] .+ (-3:0.05:3)
yrange = p[2] .+ (-3:0.05:3)
zrange = p[3] .+ (-3:0.05:3)
densities = grid_density(bluesoap_nerf, xrange, yrange, zrange)
for z in zrange
    densities = slice_density(bluesoap_nerf, xrange, yrange, z)

    plt = heatmap(1:size(densities,1),
        1:size(densities,2), densities,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x values", ylabel="y values",
        title="My title")
    display(plt)
end


points,faces = isosurface(A, MarchingCubes(iso=1))


xrange = -0.5:0.005:0.5
yrange = -0.5:0.005:0.5
zrange = -0.5:0.005:0.5
densities = grid_density(bluesoap_nerf, xrange, yrange, zrange)
maximum(densities)

xrange = -0.5:0.005:0.5
yrange = -0.5:0.005:0.5
zrange = -0.5:0.005:0.5
densities = grid_density(bluesoap_nerf, xrange, yrange, zrange)
maximum(densities)
