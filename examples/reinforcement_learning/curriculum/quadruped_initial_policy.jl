using Dojo
using IterativeLQR
using LinearAlgebra

vis = Visualizer()
open(vis)

################################################################################
# Utils
################################################################################

function load_trajectory()
	mechanism = get_quadruped()
	file = jldopen(joinpath(@__DIR__, "../../trajectory_optimization/curriculum/quadruped_4_steps_spring.jld2"))
	x_sol = file["x_sol"]
	u_sol = file["u_sol"]
	z_sol = [minimal_to_maximal(mechanism, x) for x in x_sol]
	return z_sol, u_sol
end


################################################################################
# Build Mechanism
################################################################################
opts_grad = SolverOptions(rtol=1e-4, btol=1e-4, verbose=false)
gravity = -9.81
timestep = 0.05
friction_coefficient = 0.8
damper = 1.5
spring = 20.0
ρ0 = 1e-2



mech = get_quadruped(timestep=timestep, damper=damper, spring=spring,
	gravity=gravity, body_contact=true,
	joint_limits=[[-0.5, -0.5, -2.5,],[ 0.5,  1.5, -1.0,]],
	contact=true,
	limits=true)
initialize!(mech, :quadruped, tran=[0,0,0.30], rot=[0,-0.1,0], v=[0.2,0,0])
function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = -0*[szeros(6); sones(nu-6)] * mechanism.timestep
    set_control!(mechanism, u)
    return
end
@elapsed storage = simulate!(mech, 2.0, ctrl!, record=true,
	opts=opts_grad)
visualize(mech, storage, vis=vis, show_contact=true)
################################################################################
# Initialize Policy
################################################################################
z_sol, u_sol = load_trajectory()
nu = control_dimension(mech)
nz = maximal_dimension(mech)
nz̄ = nz + 1
nzg = maximal_dimension(mech, attjac=true)
N = length(u_sol)
n = 8
NA = Int(N/n)
A = [zeros(nu,nz̄) for i = 1:NA]

function flatten_policy(A)
	nu, nz̄ = size(A[1])
	a = vcat(reshape.(A, (nu*nz̄))...)
	return a
end

function expand_policy(a; NA::Int=NA, nz̄::Int=nz+1, nu::Int=nu)
	A = [reshape(a[(i-1)nu*nz̄+1:i*nu*nz̄], (nu,nz̄)) for i=1:NA]
	return A
end

function linear_policy(z::AbstractVector, i::Int, A::AbstractVector,
	 	u_sol::AbstractVector; N::Int=N, n::Int=n)
	z̄ = [1; z]
	i1 = (i-1) % N + 1
	i2 = Int(ceil(((i-1) % N + 1) / n))
	u = u_sol[i1] + A[i2] * z̄
	return u
end

A0 = [rand(nu,nz+1) for i = 1:Int(N/n)]
a0 = flatten_policy(A0)
A00 = expand_policy(a0)
norm(A00 .- A0, Inf)


################################################################################
# Rollout policy
################################################################################
function ctrl!(mechanism, i; offset::Int=0)
	z = get_maximal_state(mechanism)
	u = linear_policy(z, i+offset, A, u_sol; N=N, n=n)
	set_control!(mechanism, u)
	return u
end
convert_frames_to_video_and_gif("quadruped_open_loop_spring")
norm(ctrl!(mech, 1) - u_sol[1], Inf)
norm(ctrl!(mech, 10) - u_sol[10], Inf)
norm(ctrl!(mech, 40) - u_sol[40], Inf)
norm(ctrl!(mech, 41) - u_sol[1], Inf)

initialize!(mech, :quadruped)
z0 = get_maximal_state(mech)
@elapsed storage = simulate!(mech, timestep*(5N+1), ctrl!, record=true,
	opts=opts_grad)
visualize(mech, storage, vis=vis, show_contact=true)
storage.x[1]

function linear_rollout(mechanism::Mechanism, A::AbstractVector; u_sol=u_sol, N::Int=N, n::Int=n)
	z = [zeros(nz) for i=1:N+1]
	dzdz = [zeros(nzg, nzg) for i=1:N]
	dzdu = [zeros(nzg, nu) for i=1:N]
	u = [zeros(nu) for i=1:N]

	initialize!(mechanism, :quadruped)
	z0 = get_maximal_state(mechanism)
	z[1] = z0
	for i = 1:N
		ctrl_offset!(m,k) = ctrl!(m,k;offset=i-1)
		simulate!(mechanism, mechanism.timestep, ctrl_offset!, record=true, opts=opts_grad)
		u[i] = linear_policy(z[i], i, A, u_sol, N=N, n=n)
		z[i+1] = get_next_state(mechanism)
		dzdz[i], dzdu[i] = get_maximal_gradients(mechanism)
		set_state!(mech, z[i+1])
	end
	return z, u, dzdz, dzdu
end

z, u, dzdz, dzdu = linear_rollout(mech, A)
visualize(mech, generate_storage(mech, z), vis=vis)

################################################################################
# Cost function and gradients
################################################################################
z = get_maximal_state(storage)[1:N+1]
u = [linear_policy(z[i], i, A, u_sol; N=N, n=n) for i = 1:N]

function tracking_cost(z, u, dzdz, dzdu, z_sol, u_sol, A; Cu=0.05, Cz=1.0)
	N = length(u_sol)
	NA = length(A)
	nz, nu = size(A[1])
	na = NA * nz * nu

	# Evaluate cost
	c = 0.0
	for i = 1:N
		c += Cz * norm(z[i+1] - z_sol[i+1])
		c += Cu * norm(u[i] - u_sol[i])
	end

	# Gradient cost
	∇c = zeros(na)




	return c, ∇c
end

function ∇tracking_cost(z, u, z_sol, u_sol, A; Cu=0.05, Cz=1.0)
	N = length(u_sol)
	NA = length(A)
	nz, nu = size(A[1])
	∇c = [zeros(nz,nu) for i=1:NA]

	dzdz = [zeros(nz,nz) for i=1:N]
	dzdu = [zeros(nz,nu) for i=1:N]
	dudA = [zeros(nz,nu) for i=1:N]
	dzdA = [zeros(nz,nu) for i=1:N]


	for i = 1:N
		dzdz[i] =
	for i = 1:N
		c += Cz * norm(z[i+1] - z_sol[i+1])
		c += Cu * norm(u[i] - u_sol[i])
	end
	return ∇c
end
A
tracking_cost(z, u, z_sol, u_sol, A)


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


function ttt(N::Int; μ0=ones(2), Σ0=Diagonal(ones(2)))
	plt = plot(aspectratio=1.0)
	X = [sqrt(Σ0) * randn(2) + μ0 for i=1:N]
	μ = sum(X) ./ N
	Xc = [x - μ for x in X]
	Σ = 1/N * sum([x*x' for x in Xc])
	Xr = [sqrt(inv(Σ)) * x for x in Xc]
	scatter!(plt, [x[1] for x in Xr], [x[2] for x in Xr], color=:blue)
	scatter!(plt, [x[1] for x in X], [x[2] for x in X], color=:red)
	display(plt)
	return Xr, μ, Σ
end

Xr, μ, Σ = ttt(400, μ0=[6,6], Σ0=Diagonal([3,3]))
Xr, μ, Σ = ttt(400, μ0=[6,6], Σ0=[5 -2; -2 5])
μ
Σ





rot = mech.joints[1].rotational
tra = mech.joints[1].translational
xa = srand(3)
xb = srand(3)
qa = UnitQuaternion(rand(4)...)
qb = UnitQuaternion(rand(4)...)
# impulse_transform_child_new(rot, xa, qa, xb, qb)
@benchmark impulse_transform_parent_new($tra, $xa, $qa, $xb, $qb)
@benchmark impulse_transform_child_new($tra, $xa, $qa, $xb, $qb)







@generated function constraint_jacobian_configuration2(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    # vec = [:(constraint_jacobian_configuration((joint.translational, joint.rotational)[$i], joint.impulses[2][λindex(joint,$i)])) for i = 1:Nc]
    vec = [:(constraint_jacobian_configuration(joint.translational, joint.impulses[2][λindex(joint,1)])) for i = 1:Nc]
    return :(cat($(vec...), dims=(1,2)))
end

function λindex2(joint::JointConstraint{T,N,Nc,RJ,TJ}, i::Int) where {T,N,Nc,RJ,TJ}
    s = 0
    for j = 1:i-1
        element = [joint.translational, joint.rotational][j]
        s += ηlength(element)
    end
    λindex2([joint.translational, joint.rotational][i], s) # to be allocation free
end


function λindex2(joint::Joint{T,Nλ,Nb,N}, s::Int) where {T,Nλ,Nb,N}
    SVector{N,Int}(s+1:s+N)
end

joint = mech.joints[2]
tr = [joint.translational, joint.rotational]
tra = joint.translational
rot = joint.rotational
i = 1
λindex2(joint, i)
Main.@code_warntype λindex2(joint, i)
Main.@code_warntype λindex2(rot, i)
Main.@code_warntype λindex2(tra, i)
@benchmark λindex2($joint, $i)
@benchmark λindex2($rot, $i)
@benchmark λindex2($tra, $i)
constraint_jacobian_configuration2(mech, joint)
@benchmark constraint_jacobian_configuration2($mech, $joint)











@benchmark fx0, fu0 = get_maximal_gradients(mech)
z = get_state(mech)
u = rand(control_dimension(mech))
fx0, fu0 = get_maximal_gradients!(mech, z, u)
fx0, fu0 = get_minimal_gradients(mech, z, u)
fx1, fu1 = get_maximal_gradients(mech)

attjac2 = cat([cat(I(6), LVᵀmat(body.state.q2[1]), I(3), dims=(1,2)) for body in mech.bodies]..., dims=(1,2))
attjac3 =
	cat(
	[cat(I(6),
	LVᵀmat(next_orientation(body.state.q2[1], body.state.ϕsol[2], mech.timestep)),
	I(3), dims=(1:2)) for body in mech.bodies]...
	, dims=(1,2))

plot(Gray.(abs.(fx0)))
plot(Gray.(abs.(attjac3' * fx1 * attjac2)))
plot(Gray.(1e6abs.(attjac3' * fx1 * attjac2 - fx0)))
norm(attjac3' * fx1 * attjac2 - fx0, Inf)


plot(Gray.(abs.(fu0)))
plot(Gray.(abs.(attjac3' * fu1)))
plot(Gray.(attjac3' * fu1 - fu0))
norm(attjac3' * fu1 - fu0, Inf)


maximal_dimension(mech, attjac=true)

z = get_state(mech)
z_next = get_next_state(mech)
x = maximal_to_minimal(mech, z)
minimal_to_maximal_jacobian(mech, x, attjac=true)
maximal_to_minimal_jacobian(mech, z_next)



full_data_matrix(mech)
nodes = [mech.joints; mech.bodies; mech.contacts]
dimrow = length.(nodes)
dimcol = data_dim.(nodes)
getfield.(nodes, :id)
indexcol = [1+sum(dimcol[1:i-1]):sum(dimcol[1:i]) for i in 1:length(dimcol)]
indexrow = [1+sum(dimrow[1:i-1]):sum(dimrow[1:i]) for i in 1:length(dimrow)]
indexcol[mech.bodies[1].id]
indexrow[mech.bodies[1].id]


datajac1 = full_matrix(mech.data_matrix, dimrow, dimcol)
full_data_matrix(mech)
plot(Gray.(datajac1))
plot(Gray.(full_data_matrix(mech)))


Main.@code_warntype displacement_jacobian_configuration(
	:parent, joint, xa, qa, xb, qb; attjac=true, vmat=true)
displacement_jacobian_configuration(
	:child, joint, xa, qa, xb, qb; attjac=true, vmat=true)
relative = :parent
@benchmark displacement_jacobian_configuration(
	$relative, $joint, $xa, $qa, $xb, $qb)
relative = :child
@benchmark displacement_jacobian_configuration(
	$relative, $joint, $xa, $qa, $xb, $qb)


a = 10
a = 10
a = 10
a = 10
a = 10