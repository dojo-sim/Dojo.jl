################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()

################################################################################
# Build Mesh
################################################################################
vis = Visualizer()
open(vis)

# NeRF density -> Mesh
rect = GeometryBasics.Rect(MeshCat.Vec(-0.75, -0.5, -0.7), MeshCat.Vec(1.5, 1, 1.4))
mesh = GeometryBasics.Mesh(x -> norm(x,0.2)-50.5, rect, MarchingCubes(), samples=(100,100,100))

# NeRF density -> density tensor -> Mesh
N = 100
xrange = range(-1.0, stop=1.0, length=N)
yrange = range(-1.0, stop=1.0, length=N)
zrange = range(-1.0, stop=1.0, length=N)
densities = grid_density(nerf_object, xrange, yrange, zrange)
vertices, faces = isosurface(densities, NaiveSurfaceNets(iso=1)) # iso=65
mesh = GeometryBasics.Mesh(densities, NaiveSurfaceNets(iso=1)) # iso=65
results_dir = joinpath(example_dir(), "results")
jldsave(joinpath(results_dir, "bunny_outer_mesh.jld2"), mesh=mesh, densities=densities,
    vertices=vertices, faces=faces)

N = 100
xrange = range(-1.0, stop=1.0, length=N)
yrange = range(-1.0, stop=1.0, length=N)
zrange = range(-1.0, stop=1.0, length=N)
densities = grid_density(nerf_object, xrange, yrange, zrange)
vertices, faces = isosurface(densities, NaiveSurfaceNets(iso=65)) # iso=65
mesh = GeometryBasics.Mesh(densities, NaiveSurfaceNets(iso=65)) # iso=65
results_dir = joinpath(example_dir(), "results")
jldsave(joinpath(results_dir, "bunny_inner_mesh.jld2"), mesh=mesh, densities=densities,
    vertices=vertices, faces=faces)

# densities = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["densities"]
# mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]


setobject!(vis[:mesh], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 1.0)))

################################################################################
# Build Convex Hull
################################################################################

using GLPK
vertices_rep = vrep(vertices)
hull_vrep = removevredundancy(vertices_rep, GLPK.Optimizer)
hull_polyhedron = polyhedron(hull_vrep)
hull_mesh = Polyhedra.Mesh(hull_polyhedron)
setobject!(vis[:hull_mesh], hull_mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))

hull_hrep = hrep(hull_polyhedron)
intersection_hrep = hull_hrep ∩ HalfSpace(SVector(0,0,1.0),0.0)
intersection_polyhedron = polyhedron(intersection_hrep)
intersection_mesh = Polyhedra.Mesh(intersection_polyhedron)

setobject!(vis[:inersection_mesh], intersection_mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))


intersection_polyhedron.vrep.points.points
for (i,p) in enumerate(intersection_polyhedron.vrep.points.points)
    obj = HyperSphere(GeometryBasics.Point(p...), 0.01)
    setobject!(vis[:vertices]["$i"], obj, MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))
end

N = 1000
samples = naive_sample(intersection_polyhedron, N)
for i = 1:N
    obj = HyperSphere(GeometryBasics.Point(samples[i,:]...), 0.01)
    setobject!(vis[:samples]["$i"], obj, MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5)))
end
