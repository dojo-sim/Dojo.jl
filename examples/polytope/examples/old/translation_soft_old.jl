################################################################################
# Translation only
################################################################################
using LinearAlgebra
using Plots
using FiniteDiff
using GLPK

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
nerf_object = py"generate_test_nerf"()

vis = Visualizer()
open(vis)

function pack_vars(p3, γ)
    return [p3; γ]
end

function unpack_vars(vars)
    p3 = vars[1:3]
    γ = vars[4:4]
    return p3, γ
end

function pack_data(p1, p2, ρ)
    return [p1; p2; ρ]
end

function unpack_data(data)
    p1 = data[1:3]
    p2 = data[4:6]
    ρ = data[7]
    return p1, p2, ρ
end

function residual(vars, data, ϕ2, ∇ϕ2; jacobian::Bool=false)
    # p1, p2, ρ = unpack_data(data)
    # p3, γ = unpack_vars(vars)
    # points2, ϕ2, ∇ϕ2 = soft_potential(p2, nerf_object)
    # points3, ϕ3, ∇ϕ3 = soft_potential(p3, nerf_object)
    # ∇ϕ2 = [0 0 1.0]
    # rdyn = m*(p3 - 2p2 + p1)/h - h*m*[0,0,g] + Vector(∇ϕ2'*γ)
    # rimp = γ - 1/ρ * [ϕ3]
    # r = [rdyn; rimp]
    # !jacobian && return r
    # Jp3 = [m/h*I(3); ∇ϕ3]
    # Jγ = [∇ϕ2 1]'
    # J = [Jp3 Jγ]
    # return r, J
    p1, p2, ρ = unpack_data(data)
    p3, γ = unpack_vars(vars)
    ∇ϕ2 = [0 0 1.0]
    rdyn = m*(p3 - 2p2 + p1)/h - h*m*[0,0,g] + (p3 - p2) * 1
    rdyn = m*(p3 - 2p2 + p1)/h - h*m*[0,0,g] - Vector(∇ϕ2'*γ) + (p3 - p2) * (1 + 100ϕ2)
    # rdyn = m*(p3 - 2p2 + p1)/h - h*m*[0,0,g] - Vector(∇ϕ2'*γ) + (p3 - p2) * ϕ2 * 100
    rimp = γ - 1/ρ * [ϕ2]
    r = [rdyn; rimp]
    !jacobian && return r
    Jp3 = [(m/h + 1 + 100ϕ2)*I(3); zeros(1,3)]
    Jγ = [∇ϕ2 1]'
    J = [Jp3 Jγ]
    return r, J
end

function get_intersection(hull_hrep, altitude)
    intersection_hrep = hull_hrep ∩ HalfSpace(SVector(0,0,1.0), -altitude)
    intersection = polyhedron(intersection_hrep)
    vrep(intersection)
    return intersection
end

function linesearch(vars, data, Δ, r, ϕ2, ∇ϕ2)
    α = 1.0
    for i = 1:10
        vars_cand = vars + α*Δ
        r_cand = residual(vars_cand, data, ϕ2, ∇ϕ2)
        (norm(r_cand, Inf) <= norm(r, Inf)) && break
        α = 1.0
    end
    return α
end

function soft_step(vars, data, ϕ2, ∇ϕ2)
    for k = 1:1
        for i = 1:10
            r, J = residual(vars, data, ϕ2, ∇ϕ2, jacobian=true)
            (norm(r, Inf) < 1e-6) && break
            # @show J
            Δ = - J \ r
            α = linesearch(vars, data, Δ, r, ϕ2, ∇ϕ2)
            vars += α * Δ
        end
        p1, p2, ρ = unpack_data(data)
        ρ /= 10.0
        data = pack_data(p1, p2, ρ)
    end
    return unpack_vars(vars)
end

function soft_simulate(vars, data, N)
    p1, p2, ρ = unpack_data(data)
    P = [p1, p2]
    Γ = []
    for t = 1:N
        points2, ϕ2, ∇ϕ2 = soft_potential(p2, nerf_object)
        p3, γ = soft_step(vars, data, ϕ2, ∇ϕ2)
        p1, p2, ρ = unpack_data(data)
        data = pack_data(p2, p3, 1.0)
        vars = pack_vars(p3, zeros(nγ))
        push!(P, p3)
        push!(Γ, γ)
    end
    return P, Γ
end

function soft_potential(x, nerf_object; N::Int=100)
    # we sample uniformly from the product of densities
        # normalized bunny density
        # soft 0-1 density for the halfspace
    altitude = x[3]
    intersection = get_intersection(hull_hrep, altitude)
    if length(intersection.vrep.points.points) > 0
        V = 0.0
        try
            V = Polyhedra.volume(intersection)
        catch e
        end
        points = naive_sample(intersection, N)
        ϕ = mean(py"density_query"(nerf_object, points)) / 1
        ∇ϕ = mean(py"density_gradient_query"(nerf_object, points), dims=1) / 1
        ϕ *= V
        ∇ϕ *= V
    else
        points = []
        ϕ = 0.0
        ∇ϕ = zeros(1,3)
    end
    return points, ϕ, ∇ϕ
end

# Nerf data
results_dir = joinpath(module_dir(), "results")
vertices = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["vertices"]
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
hull_vrep = removevredundancy(vrep(vertices), GLPK.Optimizer)
hull = polyhedron(hull_vrep)
hull_hrep = hrep(hull)

# problem parameters
m = 1.0
g = -9.81
h = 0.01

p10 = [0,0,1.0]
p20 = [0,0,1.0]
ρ0 = 1.0
data0 = pack_data(p10, p20, ρ0)

p30 = [0,0,1.0]
γ0 = [0.0]
nγ = 1
vars0 = pack_vars(p30, γ0)

p0, ϕ0, ∇ϕ0 = soft_potential([0,0,0.0], nerf_object)
H = 200
residual(vars0, data0, ϕ0, ∇ϕ0)
soft_step(vars0, data0, ϕ0, ∇ϕ0)
P0, Γ0 = soft_simulate(vars0, data0, H)


setobject!(vis[:animated], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))

animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    P0[i]
    atframe(animation, i) do
        settransform!(vis[:animated], MeshCat.Translation(P0[i]))
    end
end
setanimation!(vis, animation)

plot([i*h for i=0:H+1], [p[3] for p in P0], linewidth=3.0, xlabel="time", ylabel="altitude")



intersection_hrep = hull_hrep ∩ HalfSpace(SVector(0,0,1.0),0.10)
intersection = polyhedron(intersection_hrep)
vrep(intersection)
intersection.vrep
alt_mesh = Polyhedra.Mesh(intersection)
setobject!(vis[:alt], alt_mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))






alts = 0:0.002:1
plot(alts, [soft_potential([0,0,alt], nerf_object, N=1000)[2] for alt in alts])
plot(alts, [soft_potential([0,0,alt], nerf_object, N=100)[3][1] for alt in alts])
plot(alts, [soft_potential([0,0,alt], nerf_object, N=100)[3][3] for alt in alts])

soft_potential([0,0,1.0], nerf_object, N=10)[2]
include("../src/utils.jl")

Dojo.convert_frames_to_video_and_gif("bunny_unstable_explicit_scheme")
using Dojo
