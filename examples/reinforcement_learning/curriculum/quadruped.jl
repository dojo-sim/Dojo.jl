using Dojo
using IterativeLQR
using LinearAlgebra

vis = Visualizer()
open(vis)

# z0 = get_maximal_state(mech)
# x0 = maximal_to_minimal(mech, z0)
# x0[1:6]
# x0[7:12]
# # x0[13:14] .= -0.50 # RLhip -0.50, +0.50
# # x0[15:16] .= 1.5 # RLthigh -0.50, +1.50
# # x0[17:18] .= -2.50 # RLcalf -2.50, -1.00
# z1 = minimal_to_maximal(mech, x0)
# x1 = maximal_to_minimal(mech, z1)
# storage = generate_storage(mech, [z1])
# visualize(mech, storage, vis=vis)



gravity = -9.81
timestep = 0.05
friction_coefficient = 0.8
damper = 5.0
spring = 0.0
ρ0 = 1e-2

mech = get_quadruped(timestep=timestep, damper=damper, spring=spring,
	gravity=gravity, body_contact=true,
	joint_limits=[[-0.5, -0.5, -2.5,],[ 0.5,  1.5, -1.0,]],
	limits=true)
initialize!(mech, :quadruped)
function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = -0*[szeros(6); sones(nu-6)] * mechanism.timestep
    set_control!(mechanism, u)
    return
end
storage = simulate!(mech, 0.20, ctrl!, record=true, opts=SolverOptions(btol=1e-4, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)

file = jldopen(joinpath(@__DIR__, "../../trajectory_optimization/curriculum/quadruped_4_steps.jld2"))
x_sol = file["x_sol"]
u_sol = file["u_sol"]

nx = minimal_dimension(mech)
nu = control_dimension(mech)
N = length(u_sol)
function eval_loss(Av, x, u; α=1.0)
	l = 0.0
	A = reshape(Av, (nu,nx))
	for i = 1:N
		l += norm(A*x[i]-u[i])
	end
	return l/N + norm(Av)
end

function grad_loss(Av, x, u; α=1.0)
	ForwardDiff.gradient(Av -> eval_loss(Av, x, u, α=α), Av)
end

function hess_loss(Av, x, u; α=1.0)
	ForwardDiff.hessian(Av -> eval_loss(Av, x, u, α=α), Av)
end


function optimize(Av, x, u; α=1.0)
	for i = 1:20
		e = eval_loss(Av, x, u, α=α)
		g = grad_loss(Av, x, u, α=α)
		(norm(g, Inf) < 1e-6) && break
		H = hess_loss(Av, x, u, α=α)
		ΔAv = - H \ g
		α = linesearch(e, ΔAv, Av, x, u, α=α)
		Av += α * ΔAv
		println("i ", i, "  e", scn(e),  "  g", scn(norm(g, Inf)), "  α ", scn(α))
	end
	return Av
end

function linesearch(e, ΔAv, Av, x, u; α=1.0)
    α = 1.0
    for i = 1:8
        e_candidate = eval_loss(Av+α*ΔAv, x, u, α=α)
        (e_candidate < e) && break
        α *= 0.5
    end
    return α
end

Av = 0.01 * rand(nu*nx)
eval_loss(Av, x_sol, u_sol, α=0.05)
grad_loss(Av, x_sol, u_sol, α=0.05)
hess_loss(Av, x_sol, u_sol, α=0.05)

Av = optimize(Av, x_sol, u_sol; α=0.05)
A = reshape(Av, (nu,nx))

function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = SVector{nu}(A * get_minimal_state(mech))
	u = u_sol[k]
    set_control!(mechanism, u)
    return
end
initialize!(mech, :quadruped)
Main.@profiler storage = simulate!(mech, 1.00, ctrl!, record=true, opts=SolverOptions(btol=1e-4, verbose=false))
@elapsed storage = simulate!(mech, 2.0, ctrl!, record=true, opts=SolverOptions(btol=1e-4, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)


A * x_sol[1] - u_sol[1]
get_minimal_state(mech)

# scale x and u using ref traj and scale u by timestep
# find best A efficiently
# rollout A

# simulate quadruped efficiently


@inline function constraint(joint::Joint{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = joint_constraint(joint, xa, qa, xb, qb, η)
    e2 = minimal_coordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[SUnitRange(1,Nb½)] - (joint.joint_limits[2] .- e2);
            s[SUnitRange(Nb½+1,Nb)] - (e2 .- joint.joint_limits[1]);
            e1;
           ]
end

joint0 = mech.joints[2]
rot0 = joint0.rotational
xa = srand(3)
xb = srand(3)
qa = UnitQuaternion(rand(4)...)
qb = UnitQuaternion(rand(4)...)
η = srand(4)
# Main.@profiler constraint(rot0, xa, qa, xb, qb, η)
# @benchmark joint_constraint($rot0, $xa, $qa, $xb, $qb, $η)
# @benchmark minimal_coordinates($rot0, $xa, $qa, $xb, $qb)
@benchmark constraint($rot0, $xa, $qa, $xb, $qb, $η)
