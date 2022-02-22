# ENV["PYCALL_JL_RUNTIME_PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"
# ENV["PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"

using BenchmarkTools
using Plots
using Pkg
using PyCall
# Pkg.build("PyCall")

osf_path = joinpath("/home/simon/research/repos/osf-pytorch")
# pushfirst!(pyimport("sys")."path", "")
pushfirst!(pyimport("sys")."path", osf_path)


################################################################################
# Query NeRf
################################################################################
@pyinclude(joinpath(osf_path, "extract_density_julia.py"))
render_kwargs_test = py"generate_test_nerf"()
xyz = rand(Float32, 100,3) ./ 5
density = py"density_query"(render_kwargs_test, xyz)
density_grad = py"density_gradient_query"(render_kwargs_test, xyz)
@belapsed density = py"density_query"(render_kwargs_test, xyz)
@belapsed density_grad = py"density_gradient_query"(render_kwargs_test, xyz)

# evaluation & gradient benchmark
Ns = [1,3,10,30,100,300,1000,3000,10000]
xyzs = [rand(Float32, N,3) for N in Ns]
eval_times = [@elapsed py"density_query"(render_kwargs_test, xyz) for xyz in xyzs]
grad_times = [@elapsed py"density_gradient_query"(render_kwargs_test, xyz) for xyz in xyzs]
plot(Ns, eval_times, axis=:log, label="evaluation", xlabel="number of samples", ylabel="time (s)",
    title="batch querying is more efficient")
plot!(Ns, grad_times, axis=:log, label="gradient", xlabel="number of samples", ylabel="time (s)",
    title="batch querying is more efficient")



################################################################################
# Build Scan
################################################################################
function z_slice(xrange, yrange, z)
    nx = length(xrange)
    ny = length(yrange)
    xv = Vector(xrange)
    yv = Vector(yrange)
    xyz = zeros(Float32, nx*ny, 3)
    for i = 1:nx
        for j = 1:ny
            xyz[(i-1)*ny+j,:] .= [xv[i], yv[j], z]
        end
    end
    return xyz
end

function density_slice(xrange, yrange, z)
    xyz = z_slice(xrange, yrange, z)
    density = py"density_query"(render_kwargs_test, xyz)
    nx = length(xrange)
    ny = length(yrange)
    density = reshape(density, nx, ny)
    return density
end

xrange = -1.0:0.01:1.0
yrange = -1.0:0.01:1.0
z = 0.1
density = density_slice(xrange, yrange, z)

heatmap(xrange, yrange, density,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x", ylabel="y",
    title="bunny")

zrange = -0.5:0.01:0.65
anim = @animate for z ∈ zrange
    xyz = z_slice(xrange, yrange, z)
    density = py"density_query"(render_kwargs_test, xyz)
    density = reshape(density, nx, ny)
    heatmap(xrange, yrange, density,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x", ylabel="y",
        title="bunny")
end

nerf_dir = joinpath(module_dir(), "examples/nerf")
# gif(anim, joinpath(nerf_dir, "bunny_scan.gif"), fps=15)


################################################################################
# Build Mesh
################################################################################
using Meshing
using MeshCat
using GeometryBasics
using LinearAlgebra
using JLD2

# vis = Visualizer()
# open(vis)

# NeRF density -> Mesh
rect = GeometryBasics.Rect(MeshCat.Vec(-0.75, -0.5, -0.7), MeshCat.Vec(1.5, 1, 1.4))
mesh = GeometryBasics.Mesh(x -> norm(x,0.2)-50.5, rect, MarchingCubes(), samples=(100,100,100))
# mesh = GeometryBasics.Mesh(
#     x -> py"density_query"(render_kwargs_test, [x[1] x[2] x[3]])[1],
#     rect, NaiveSurfaceNets(iso=65), samples=(25,25,25))
# jldsave(joinpath(nerf_dir, "bunny_mesh_.obj"), mesh=mesh, D=D)

mesh = jldopen(joinpath(nerf_dir, "bunny_mesh.obj"))["mesh"]
setobject!(vis, mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))

# NeRF density -> density tensor -> Nesh
N = 100
xrange = range(-1.0, stop=1.0, length=N)
yrange = range(-1.0, stop=1.0, length=N)
zrange = range(-1.0, stop=1.0, length=N)
# D = zeros(N,N,N)
# for i = 1:N
#     D[:,:,i] = density_slice(xrange, yrange, zrange[i])
# end

mesh_old = jldopen(joinpath(nerf_dir, "bunny_mesh_low_res.obj"))["mesh"]
mesh_new = jldopen(joinpath(nerf_dir, "bunny_mesh.obj"))["mesh"]
D = jldopen(joinpath(nerf_dir, "bunny_mesh.obj"))["D"]
points,faces = isosurface(D, NaiveSurfaceNets(iso=65))
mesh = GeometryBasics.Mesh(D, NaiveSurfaceNets(iso=65))
setobject!(vis[:old], mesh_old,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))
setobject!(vis[:new], mesh_new,
    MeshPhongMaterial(color=RGBA{Float32}(0.5, 0.5, 1, 0.5)))
