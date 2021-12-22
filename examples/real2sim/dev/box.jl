# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "real2sim", "utils.jl"))

mech = getmechanism(:box, Δt=0.05, g=-9.81, cf=0.2, radius=0.00, side=0.50);
initialize!(mech, :box, x=[0,-1,1.], v=[0,2,1.], ω=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact = false)


body1 = mech.bodies.values[1]
ineqc1 = mech.ineqconstraints.values[1]
bnd = ineqc1.constraints[1]

Δt = mech.Δt
q2 = body1.state.q2[1]
ϕ25 = body1.state.ϕsol[2]
x3, q3 = posargs3(body1.state, Δt)
γ = ineqc1.γsol[2]

function d(vars)
	x = vars[1:3]
	q = UnitQuaternion(vars[4:7]..., false)
	return ∂g∂ʳpos(bnd, x, q, nothing)' * γ
end

M = ∂integration(q2, ϕ25, Δt)
D0 = FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M

X = [bnd.ainv3;
     szeros(1,3);
     bnd.Bx]
∇ = -X * VLmat(q) * RᵀVᵀmat(q) * skew(bnd.p - vrotate(bnd.offset, inv(q3)))

D1 = [szeros(Float64,6,3) [szeros(Float64,3,3); ∇]]


q = UnitQuaternion(rand(4)...)
p = rand(3)
offset = rand(3)
λ = rand(3)


Q = - X * q * skew(p - vrotate(offset, inv(q)))
res0 = Q'*γ

Q = - X * VLmat(q) * RᵀVᵀmat(q) * skew(p - vrotate(offset, inv(q)))
res1 = Q'*γ

res2 = skew(p - vrotate(offset, inv(q))) * (VLmat(q) * RᵀVᵀmat(q))' * X' * γ
res3 = skew(p - vrotate(offset, inv(q))) * VRmat(q) * LᵀVᵀmat(q) * X' * γ

norm(res0 - res1, Inf)
norm(res0 - res2, Inf)
norm(res0 - res3, Inf)


res_fct(q, p, offset, λ) = skew(p - vrotate(offset, inv(q))) * VRmat(q) * LᵀVᵀmat(q) * λ
# res_fct(q, p, offset, λ) = VRmat(q) * LᵀVᵀmat(q) * λ
# res_fct(q, p, offset, λ) = skew(p - vrotate(offset, inv(q))) * λ
Dres0 = FiniteDiff.finite_difference_jacobian(q -> res_fct(UnitQuaternion(q..., false), p, offset, λ), vector(q))

# Dres1 = VRmat(q) * ∂qLᵀVᵀmat(λ) + ∂qVRmat(LᵀVᵀmat(q) * λ)
# Dres1 = ∂pskew(λ) * -∂vrotate∂q(offset, inv(q)) * Tmat()
Dres1 = ∂pskew(VRmat(q) * LᵀVᵀmat(q) * λ) * -∂vrotate∂q(offset, inv(q)) * Tmat()
Dres1 += skew(p - vrotate(offset, inv(q))) * ∂qVRmat(LᵀVᵀmat(q) * λ)
Dres1 += skew(p - vrotate(offset, inv(q))) * VRmat(q) * ∂qLᵀVᵀmat(λ)


norm(Dres0 - Dres1, Inf)



Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat()

Vmat() * Rmat(q) * Lmat(q)' * Vᵀmat()

Qfct(q,p,offset) = skew(p - vrotate(offset, inv(q))) * q * X' * γ


using Symbolics
@variables p[1:3]
@variables λ[1:3]
@variables offset[1:3]
@variables k[1:3]
@variables l[1:4]
@variables q[1:4]

qq = UnitQuaternion(q..., false)
myexpr = skew(p - vrotate(offset, inv(qq))) * VRmat(qq) * LᵀVᵀmat(qq) * λ
myjacexpr = Symbolics.jacobian(myexpr, q)


skew(k) * λ
skew(λ)
RᵀVᵀmat(qq) * k
LᵀVᵀmat(qq) * k
VLmat(qq) * l
VRmat(qq) * l
RᵀVᵀexpr = Symbolics.jacobian(RᵀVᵀmat(qq) * k, q)
LᵀVᵀexpr = Symbolics.jacobian(LᵀVᵀmat(qq) * k, q)
VLexpr = Symbolics.jacobian(VLmat(qq) * l, q)
VRexpr = Symbolics.jacobian(VRmat(qq) * l, q)
skewexpr = Symbolics.jacobian(skew(k) * λ, k)







q = UnitQuaternion(rand(4)...)
r = rand(3)
R0 = hmat()' * lmult(q) * rmult(q)' * hmat() * r
R1 = Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * r
norm(R0 - R1, Inf)

hmat()
Vmat()
using Symbolics



################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [-2ones(3), [2,2,-1.]],
				   :ωlims => [-6ones(3), 6ones(3)])
mech_kwargs = Dict(:cf => 0.1, :radius => 0.00, :side => 0.5)
generate_dataset(:box, H=0.40, N=25,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	sleep_ratio=1,
	show_contact=false)


################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:box; N = 25, mech_kwargs...)

data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
global const ROTATE = Ref{Float64}(0.0)
@profiler loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[loss(:box, pairs0, data0 + [i;zeros(39)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[loss(:box, pairs0, data0 + [0;i;zeros(38)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################
using Optim

function termination_callback(opt; value_tol=1e-4)
	return opt.value < value_tol
end

solver = LBFGS(;m=100,
        alphaguess=Optim.LineSearches.InitialStatic(),
        linesearch=Optim.LineSearches.HagerZhang(),
        P = 3e1*I(25),
		)

lower = [0.00,
	+0.05, +0.05, -1.00,
	+0.05, -1.00, -1.00,
	-1.00, +0.05, -1.00,
	-1.00, -1.00, -1.00,
	+0.05, +0.05, +0.05,
	+0.05, -1.00, +0.05,
	-1.00, +0.05, +0.05,
	-1.00, -1.00, +0.05]
upper = [0.80,
	+1.00, +1.00, -0.05,
	+1.00, -0.05, -0.05,
	-0.05, +1.00, -0.05,
	-0.05, -0.05, -0.05,
	+1.00, +1.00, +1.00,
	+1.00, -0.05, +1.00,
	-0.05, +1.00, +1.00,
	-0.05, -0.05, +1.00]

function d2data(d)
	cf = d[1]
	data = [cf; 0; +d[2:4];
			cf; 0; +d[5:7];
			cf; 0; +d[8:10];
			cf; 0; +d[11:13];
			cf; 0; +d[14:16];
			cf; 0; +d[17:19];
			cf; 0; +d[20:22];
			cf; 0; +d[23:25];
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(25))

function fg!(F,G,d)
	# do common computations here
	l, ∇ = loss(:box, pairs0, d2data(d), n_sample=50)
	if G != nothing
		G .= ∇d2data' * ∇
	end
	if F != nothing
		value = l
	    return value
	end
end

d0 = [0.40,
	+0.50, +0.50, -0.50,
	+0.50, -0.50, -0.50,
	-0.50, +0.50, -0.50,
	-0.50, -0.50, -0.50,
	+0.50, +0.50, +0.50,
	+0.50, -0.50, +0.50,
	-0.50, +0.50, +0.50,
	-0.50, -0.50, +0.50]

ROTATE[] = 0.0
optimize(Optim.only_fg!(fg!), lower, upper, d0, Fminbox(solver),
	Optim.Options(
		callback = termination_callback,
		# allow_f_increases = true,
		show_trace = true),
	; inplace = false,
	)



# function callback(optim_state)
# 	# @show fieldnames(typeof(optim_state.metadata))
# 	# Main.SEED += 1/40
# 	return false
# end

# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS


################################################################################
# Visualization
################################################################################
