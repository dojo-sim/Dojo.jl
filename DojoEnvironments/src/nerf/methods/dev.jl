


v = Vector(-1:0.01:1)
plot([coulomb_direction(x, 1e1, 1e-8)[1] for x in v])
coulomb_direction(-10:0.01:10, 1e3)
# Set data
Nb = length(mech.bodies)
data = Dojo.get_data(mech)
Dojo.set_data!(mech, data)
sol = Dojo.get_solution(mech)
attjac = Dojo.attitude_jacobian(data, Nb)

# IFT
solmat = -Dojo.full_matrix(mech.system)
# finite diff
fd_solmat = finite_difference_solution_matrix(mech, data, sol,
    δ=1.0e-8)

plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))
plot(Gray.(100000abs.(solmat - fd_solmat)))
norm(solmat - fd_solmat, Inf)
norm(solmat[51:53,1:12] - fd_solmat[51:53,1:12], Inf)

solmat[51:53,1:12]
fd_solmat[51:53,1:12]

mech.joints
mech.bodies
mech.contacts






using FiniteDiff
using Plots

timestep = 0.01
model = mech.contacts[1].model
model.collision.collider.options.impact_spring = 1e4
model.collision.collider.options.impact_damper = 2e4
model.collision.collider.options.sliding_drag = 3e0
model.collision.collider.options.sliding_friction = 4e0
model.collision.collider.options.rolling_drag = 1e0
model.collision.collider.options.rolling_friction = 1e0
model.collision.collider.options.coulomb_smoothing = 3e0
model.collision.collider.options.coulomb_regularizer = 1e-0

vp = srand(3)
xp = srand(3)
qp = rand(4)
qp = Quaternion(qp ./ norm(qp)..., true)
ϕp = srand(3)

xc = srand(3)
vc = srand(3)
qc = rand(4)
qc = Quaternion(qc ./ norm(qc)..., true)
ϕc = srand(3)

x2p = next_position(xp, -vp, timestep)
q2p = next_orientation(qp, -ϕp, timestep)
x2c = next_position(xc, -vc, timestep)
q2c = next_orientation(qc, -ϕc, timestep)

function local_constraint(model, x2p, vp, q2p, ϕp, x2c, vc, q2c, ϕc, timestep)
    xp = next_position(x2p, vp, timestep)
    qp = next_orientation(q2p, ϕp, timestep)
    xc = next_position(x2c, vc, timestep)
    qc = next_orientation(q2c, ϕc, timestep)
    return constraint(model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
end

J0 = constraint_jacobian_velocity(:parent, model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> local_constraint(model,
    x2p, vϕ[SUnitRange(1,3)], q2p, vϕ[SUnitRange(4,6)], x2c, vc, q2c, ϕc, timestep),
    [vp; ϕp])
norm(J0 - J1, Inf)
plot(Gray.(100abs.(Matrix(J0 - J1))))


v = srand(3)
J0 = 1/norm(v) * Diagonal(sones(3)) - v*v' ./ norm(v)^3
J1 = FiniteDiff.finite_difference_jacobian(v -> v/norm(v), v)
J0 - J1






















deps_folder = joinpath(module_dir(), "environments/bunny/deps")
inner_mesh_path = joinpath(deps_folder, "bunny_inner_mesh.obj")
outer_mesh_path = joinpath(deps_folder, "bunny_outer_mesh.obj")
collider = jldopen(joinpath(deps_folder, "soft_collider.jld2"))["soft"]
soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path)
[soft_body,]



mech.bodies[1].shape
Body(1.0, szeros(3,3), shape=mech.bodies[1].shape)
collider = mech.contacts[1].model.collider
bodies = [SoftBody(collider, name=:bunny, color=RGBA(1,1,1,0.1))]

# shape0 = bodies[1].shape
shape1 = Mesh(joinpath(deps_folder, "bunny_inner_mesh.obj"),
    position_offset = -collider.center_of_mass,
    color=RGBA(0.2,0.2,0.2,1.0))
shape2 = Mesh(joinpath(deps_folder, "bunny_outer_mesh.obj"),
    position_offset = -collider.center_of_mass,
    color=RGBA(0.9,0.9,0.9,0.3))
shape_vec = Vector{Shape{T}}([shape1, shape2])#, shape0]
shapes = Shapes(shape_vec)
bodies[1].shape = shapes



using Plots
using Dojo

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
light_blue = RGBA(0.4,0.4,1.0,0.8)

f(x) = x
ϕ(x) = x
L(x, ρ) = f(x) - ρ * log.(ϕ(x))
X = 0:0.0001:1
Y = -1:0.1:1
plt = plot(ylims=(-0.2,1.0), xlims=(-0.2,1.0), yticks=[0,1], xticks=[0,1], legend=:bottomright, size=(300,300))
plot!(plt, Y, f.(Y), linewidth=6.0, color=:black, label="f(x)")
plot!(rectangle(-0.2,2,0,-1), opacity=.6, color=:red, label="ϕ(x) < 0")
anim = @animate for i = 1:10
    plot!(plt, X, [1; L.(X[2:end], 10*0.5^i)], linewidth=5.0, color=light_blue, label=false)# label="ρ = 3e-1")
end
gif(anim, fps=10, "/home/simon/Downloads/anim.gif")

PHI = Vector(1:-0.005:0)
plt = plot(ylims=(-0.01,1.0), xlims=(-0.01,1.0), yticks=[0,1], xticks=[0,1], legend=:topright, size=(300,300), )
anim = @animate for i = 1:440
    # for (j,κ) ∈ enumerate([1e-1, 5e-2, 1e-2])
    for (j,κ) ∈ enumerate([5e-2,])
        PHI_i = PHI[1:min(i,length(PHI))]
        (i == 1) && plot!(plt, PHI_i, κ ./ PHI_i, color=light_blue, linewidth=1+1.5j, label="κ = $κ")
        plot!(plt, PHI_i, κ ./ PHI_i, color=light_blue, linewidth=1+1.5j, label=nothing)
    end
end

gif(anim, fps=40, "/home/simon/Downloads/force_solo.gif")

plot!(plt, PHI, 1e-2 ./ PHI, color=light_blue, linewidth=5.0, label="κ = 1e-2")
plot!(plt, PHI, 1e-3 ./ PHI, color=light_blue, linewidth=5.0, label="κ = 1e-3")






J0, J1 = test_solmat(:bunny_triumvirate)




function test_solmat(model;
    ϵ=1.0e-6,
    tsim=0.1,
    ctrl=(m, k)->nothing,
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    verbose=false,
    T=Float64,
    kwargs...)

    # mechanism
    mechanism = get_mechanism(model,
        timestep=timestep,
        gravity=gravity;
        kwargs...)
    initialize!(mechanism, model)

    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true,
        verbose=verbose,
        opts=SolverOptions(rtol=ϵ, btol=ϵ))

    # Set data
    Nb = length(mechanism.bodies)
    data = Dojo.get_data(mechanism)
    Dojo.set_data!(mechanism, data)
    sol = Dojo.get_solution(mechanism)
    attjac = Dojo.attitude_jacobian(data, Nb)

    # IFT
    solmat = Dojo.full_matrix(mechanism.system)
    # finite diff
    fd_solmat = finite_difference_solution_matrix(mechanism, data, sol,
        δ=1.0e-5,
        verbose=verbose)
    return fd_solmat, solmat
end

function finite_difference_solution_matrix(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector;
    δ=1.0e-8,
    verbose=false)

    nsol = length(sol)
    jac = zeros(nsol, nsol)

    Dojo.set_data!(mechanism, data)
    Dojo.set_solution!(mechanism, sol)

    for i = 1:nsol
        verbose && println("$i / $nsol")
        solp = deepcopy(sol)
        solm = deepcopy(sol)
        solp[i] += δ
        solm[i] -= δ
        rp = Dojo.evaluate_residual!(deepcopy(mechanism), data, solp)
        rm = Dojo.evaluate_residual!(deepcopy(mechanism), data, solm)
        jac[:, i] = (rp - rm) / (2δ)
    end
    return jac
end

function control!(mechanism, k;
    u=0.1)
    for joint in mechanism.joints
        nu = Dojo.input_dimension(joint,
            ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        Dojo.set_input!(joint, su)
    end
end
