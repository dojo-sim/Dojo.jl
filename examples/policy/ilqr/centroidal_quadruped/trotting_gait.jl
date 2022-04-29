# ## model
include("trajopt_model.jl")

# ## horizon

freq = 100
h = 1.0 / freq
T = Int(floor(0.5 * 0.65 / h)) + 1
Tm = 17 # T / 2

# ## centroidal_quadruped
s = get_simulation("centroidal_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env

nx = 2 * model.nq
nc = 4 #model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 5

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(model, env, [h], y, x, u, w), nx + nθ + nx + model.nu, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(model, env, [h], y, x, u, w), nx + nθ + nx + model.nu, nx + nθ + nx + model.nu, nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
mode = :left
body_height = 0.3
foot_x = 0.17
foot_y = 0.15
foot_height = 0.08

q1 = zeros(model.nq)

# body position
q1[1:3] = [0.0; 0.0; body_height]
q1[4:6] = [0.0; 0.0; 0.0]

# foot1
q1[7:9] = [foot_x; foot_y; 0]

# foot2
q1[10:12] = [foot_x; -foot_y; 0]

# foot3
q1[13:15] = [-foot_x; foot_y; 0]

# foot4
q1[16:18] = [-foot_x; -foot_y; 0]

qM = deepcopy(q1)

if mode == :left
    qM[6 + 3] += foot_height # front left
    qM[15 + 3] += foot_height # back right
elseif mode == :right
    qM[9 + 3] += foot_height # front right
    qM[12 + 3] += foot_height # back left
elseif mode == :none
end

qT = copy(q1)

# visualize!(vis, model, [qM], Δt=h);

q_ref = [linear_interpolation(q1, qM, Tm)..., linear_interpolation(qM, qT, Tm)...]
x_ref = [[q_ref[t]; q_ref[t+1]] for t = 1:T]
x1 = x_ref[1]
xM = x_ref[Tm]
xT = x_ref[T]

visualize!(vis, model, q_ref, Δt=h);

# ## objective
obj = DTO.Cost{Float64}[]
for t = 1:T
    if t == T
        function objT(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-1 * dot(v, v)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            return J
        end
        push!(obj, DTO.Cost(objT, nx + nθ + nx + model.nu, 0))
    elseif t == 1
        function obj1(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            J += 0.5 * transpose(u[1:model.nu]) * Diagonal(1.0e-3 * ones(model.nu)) * u[1:model.nu]
            J += 1000.0 * u[end] # slack
            return J
        end
        push!(obj, DTO.Cost(obj1, nx, nu))
    else
        function objt(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            u_previous = x[nx + 5 + nx .+ (1:model.nu)]
            u_control = u[1:model.nu]
            w = (u_control - u_previous) ./ h
            J += 0.5 * 1.0e-2 * dot(w, w)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            J += 0.5 * transpose(u[1:model.nu]) * Diagonal(1.0e-3 * ones(model.nu)) * u[1:model.nu]
            J += 1000.0 * u[end] # slack
            return J
        end
        push!(obj, DTO.Cost(objt, nx + nθ + nx + model.nu, nu))
    end
end

# ## constraints
# initial condition
xl1 = [q1; q1]
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)]
xut = [Inf * ones(nx); Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

# final condition
xlT = [qT; qT; -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)]
xuT = [qT; qT; Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
bndt = DTO.Bound(nx + nθ + nx + model.nu, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
bndT = DTO.Bound(nx + nθ + nx + model.nu, 0, xl=xlT, xu=xuT)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

cons = DTO.Constraint{Float64}[]
for t = 1:T
    if t == 1
        function constraints_1(x, u, w)
            [
            # equality (16)
            contact_constraints_equality(model, env, h, x, u, w);
            # inequality (28)
            contact_constraints_inequality_1(model, env, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28))))
    elseif t == T
        function constraints_T(x, u, w)
            [
            # inequality (8)
            contact_constraints_inequality_T(model, env, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_T, nx + nθ + nx + model.nu, nu, idx_ineq=collect(0 .+ (1:8))))
    else
        function constraints_t(x, u, w)
            [
            # equality (16)
            contact_constraints_equality(model, env, h, x, u, w);
            # inequality (32)
            contact_constraints_inequality_t(model, env, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_t, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32))) )
    end
end

# ## problem
tolerance = 1.0e-3
p = DTO.solver(dyn, obj, cons, bnds,
    options=DTO.Options(
        max_iter=2000,
        tol=tolerance,
        constr_viol_tol=tolerance,
        ))

# ## initialize
x_interpolation = [x_ref[1], [[x_ref[t]; zeros(nθ); zeros(nx); zeros(model.nu)] for t = 2:T]...]
u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(p, x_interpolation)
DTO.initialize_controls!(p, u_guess)

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show x_sol[1]
@show x_sol[T]
sum([u[end] for u in u_sol[1:end-1]])

# ## visualize
vis = Visualizer()
render(vis)
visualize!(vis, model, x_sol, Δt=h);

q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
u_opt = [u[1:model.nu] for u in u_sol]
λ_opt = [u[model.nu .+ (1:4)] for u in u_sol]
b_opt = [u[model.nu + 4 .+ (1:16)] for u in u_sol]

# # mirror
# q_mirror = Vector{Float64}[]
# u_mirror = Vector{Float64}[]
# reflection = Diagonal([1.0; -1.0; 1.0])

# for (t, q) in enumerate(q_opt)
#     qfl = q[6 .+ (1:3)]  # front left
#     qfr = q[9 .+ (1:3)]  # front right
#     qbl = q[12 .+ (1:3)] # back left
#     qbr = q[15 .+ (1:3)] # back right

#     push!(q_mirror, [q[1:6]; reflection * qfr; reflection * qfl; reflection * qbr; reflection * qbl])

#     if t < T
#         # set left using right
#         ufl = u_opt[t][0 .+ (1:3)] # front left
#         ufr = u_opt[t][3 .+ (1:3)] # front right
#         ubl = u_opt[t][6 .+ (1:3)] # back left
#         ubr = u_opt[t][9 .+ (1:3)] # back right
#         if t < T
#             push!(u_mirror, [reflection * ufr; reflection * ufl; reflection * ubr; reflection * ubl])
#         end
#     end
# end

# using Plots
# plot(hcat(q_opt..., q_mirror...)', xlabel="time", ylabel="configuration", label="")
# plot(hcat(u_opt..., u_mirror...)', xlabel="time", ylabel="control", label="")

# RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, [q_opt..., q_mirror...], Δt=h);

# function mirror_gait(q, u, γ, b, ψ, η, T)
# 	qm = [deepcopy(q)...]
# 	um = [deepcopy(u)...]
# 	γm = [deepcopy(γ)...]
# 	bm = [deepcopy(b)...]
# 	ψm = [deepcopy(ψ)...]
# 	ηm = [deepcopy(η)...]

#     # stride = zero(qm[1])
# 	@show strd = q[T+1][1] - q[2][1]

# 	p10 = q[2][6 .+ (1:3)]
# 	p20 = q[2][9 .+ (1:3)]
# 	p30 = q[2][12 .+ (1:3)]
# 	p40 = q[2][15 .+ (1:3)]

# 	for t = 1:T-1
# 		# configuration
# 		pt = q[t+2][1:3]
# 		a = q[t+2][3 .+ (1:3)]
# 		p1 = q[t+2][6 .+ (1:3)]
# 		# p2 = q[t+2][9 .+ (1:3)]
# 		# p3 = q[t+2][12 .+ (1:3)]
# 		p4 = q[t+2][15 .+ (1:3)]

# 		p2_diff = q[t+2][9 .+ (1:3)] - p20
# 		p3_diff = q[t+2][12 .+ (1:3)] - p30

# 		push!(qm, [pt + [0.5 * strd; 0.0; 0.0]; a;
# 			p10[1] + p2_diff[1]; q[t+2][8]; p10[3] + p2_diff[3]
# 			q[end][9 .+ (1:3)];
# 			q[end][12 .+ (1:3)];
# 			p40[1] + p3_diff[1]; q[t+2][17]; p40[3] + p3_diff[3];
# 			])
# 		# push!(qm, q[t+2])

# 		# control
# 		u1 = u[t][1:3]
# 		u2 = u[t][3 .+ (1:3)]
# 		u3 = u[t][6 .+ (1:3)]
# 		u4 = u[t][9 .+ (1:3)]
# 		push!(um, [u2; u1; u4; u3])

# 		# impact
# 		γ1 = γ[t][1]
# 		γ2 = γ[t][2]
# 		γ3 = γ[t][3]
# 		γ4 = γ[t][4]
# 		push!(γm, [γ2; γ1; γ4; γ3])
# 		# push!(γm, γ[t])

# 		# friction
# 		b1 = b[t][1:4]
# 		b2 = b[t][4 .+ (1:4)]
# 		b3 = b[t][8 .+ (1:4)]
# 		b4 = b[t][12 .+ (1:4)]
# 		push!(bm, [b2; b1; b4; b3])
# 		# push!(bm, b[t])

# 		# dual
# 		ψ1 = ψ[t][1]
# 		ψ2 = ψ[t][2]
# 		ψ3 = ψ[t][3]
# 		ψ4 = ψ[t][4]
# 		push!(ψm, [ψ2; ψ1; ψ4; ψ3])
# 		# push!(ψm, ψ[t])

# 		# dual
# 		η1 = η[t][1:4]
# 		η2 = η[t][4 .+ (1:4)]
# 		η3 = η[t][8 .+ (1:4)]
# 		η4 = η[t][12 .+ (1:4)]
# 		push!(ηm, [η2; η1; η4; η3])
# 		# push!(ηm, η[t])
# 	end

# 	return qm, um, γm, bm, ψm, ηm
# end

q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
u_opt = [u[1:model.nu] for u in u_sol]
γ_opt = [u[model.nu .+ (1:4)] for u in u_sol]
b_opt = [u[model.nu + 4 .+ (1:16)] for u in u_sol]
ψ_opt = [u[model.nu + 4 + 16 .+ (1:4)] for u in u_sol]
η_opt = [u[model.nu + 4 + 16 + 4 .+ (1:16)] for u in u_sol]

if mode == :right
    qr = q_opt
    vr = v_opt
    ur = u_opt
    γr = γ_opt
    br = b_opt
    ψr = ψ_opt
    ηr = η_opt

    using JLD2
    @save joinpath(@__DIR__, "inplace_trot_right.jld2") qr ur γr br ψr ηr
    @load joinpath(@__DIR__, "inplace_trot_right.jld2") qr ur γr br ψr ηr
end

if mode == :left
    ql = q_opt
    vl = v_opt
    ul = u_opt
    γl = γ_opt
    bl = b_opt
    ψl = ψ_opt
    ηl = η_opt

    using JLD2
    @save joinpath(@__DIR__, "inplace_trot_left.jld2") ql ul γl bl ψl ηl
    @load joinpath(@__DIR__, "inplace_trot_left.jld2") ql ul γl bl ψl ηl
end

if mode == :none
    ql = q_opt
    vl = v_opt
    ul = u_opt
    γl = γ_opt
    bl = b_opt
    ψl = ψ_opt
    ηl = η_opt

    using JLD2
    @save joinpath(@__DIR__, "inplace_trot_none.jld2") ql ul γl bl ψl ηl
    @load joinpath(@__DIR__, "inplace_trot_none.jld2") ql ul γl bl ψl ηl
end

@load joinpath(@__DIR__, "inplace_trot_left.jld2") ql ul γl bl ψl ηl

qm = [qr[2:end]..., ql[2:end]...]
um = [ur..., ul...]
γm = [γr..., γl...]
bm = [br..., bl...]
ψm = [ψr..., ψl...]
ηm = [ηr..., ηl...]
μm = model.μ_world
hm = h
timesteps = range(0.0, stop=(h * (length(qm) - 2)), length=(length(qm) - 2))
plot(timesteps, hcat(qm[2:end-1]...)', labels="")
plot(timesteps, hcat(um...)', labels="")
plot(timesteps, hcat(γm...)', labels="")
plot(timesteps, hcat(bm...)', labels="")
plot(timesteps, hcat(ψm...)', labels="")
plot(timesteps, hcat(ηm...)', labels="")

using JLD2
@save joinpath(@__DIR__, "inplace_trot.jld2") qm um γm bm ψm ηm μm hm
@load joinpath(@__DIR__, "inplace_trot.jld2") qm um γm bm ψm ηm μm hm
