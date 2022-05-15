using OSFLoader
using GeometryBasics
using LinearAlgebra
using JLD2
using MeshIO
using Meshing
using MeshCat
using Quaternions

results_dir = joinpath(OSFLoader.osf_loader_dir(), "assets/collider")
num_particle = 5000

################################################################################
# bunny
################################################################################
bunny_nerf = get_nerf_object(filename="bunny_trans")
bunny_normalizer = DensityFieldNormalizer(nerf=:bunny)
mass, inertia, center_of_mass = inertia_properties(bunny_nerf; normalizer=bunny_normalizer)
particles, densities, density_gradients = sample_soft(bunny_nerf, num_particle,
    normalizer=bunny_normalizer, particle_noise=0.005)

JLD2.jldsave(joinpath(results_dir, "bunny.jld2"),
    mass=mass,
    inertia=inertia,
    center_of_mass=center_of_mass,
    particles=particles,
    densities=densities,
    density_gradients=density_gradients,
    )


################################################################################
# bluesoap
################################################################################
bluesoap_nerf = get_nerf_object(filename="bluesoap")
bluesoap_normalizer = DensityFieldNormalizer(nerf=:bluesoap)
mass, inertia, center_of_mass = inertia_properties(bluesoap_nerf; normalizer=bluesoap_normalizer)
particles, densities, density_gradients = sample_soft(bluesoap_nerf, num_particle,
    normalizer=bluesoap_normalizer, particle_noise=0.005)

JLD2.jldsave(joinpath(results_dir, "bluesoap.jld2"),
    mass=mass,
    inertia=inertia,
    center_of_mass=center_of_mass,
    particles=particles,
    densities=densities,
    density_gradients=density_gradients,
    )


################################################################################
# halfsoap
################################################################################
halfsoap_nerf = get_nerf_object(filename="halfsoap")
halfsoap_normalizer = DensityFieldNormalizer(nerf=:halfsoap)
mass, inertia, center_of_mass = inertia_properties(halfsoap_nerf; normalizer=halfsoap_normalizer)
particles, densities, density_gradients = sample_soft(halfsoap_nerf, num_particle,
    normalizer=halfsoap_normalizer, particle_noise=0.005)

JLD2.jldsave(joinpath(results_dir, "halfsoap.jld2"),
    mass=mass,
    inertia=inertia,
    center_of_mass=center_of_mass,
    particles=particles,
    densities=densities,
    density_gradients=density_gradients,
    )
