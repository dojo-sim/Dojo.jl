# skew-symmetric matrix
function hat(x)
    return [0 -x[3] x[2];
            x[3] 0 -x[1];
           -x[2] x[1] 0]
end

function inverse_quaternion(q) 
	T = Diagonal([1.0; -1; -1; -1])
	return T * q 
end

# left quaternion multiply as matrix
function L_mult(x)
    [x[1] -transpose(x[2:4]); 
     x[2:4] x[1] * I(3) + hat(x[2:4])]
end

# right quaternion multiply as matrix
function R_mult(x)
    [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - hat(x[2:4])]
end

# https://roboticexplorationlab.org/papers/planning_with_attitude.pdf
function attitude_jacobian(x)
    H = [zeros(1, 3); I(3)]
    return L_mult(x) * H
end

# rotation matrix
function rotation_matrix(q) 
    H = [zeros(1, 3); I(3)]
    transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
end

# right discrete Legendre transform for a free rigid body
function DLT1_linear(h, m, x1, x2)
	-1.0 / h * m * (x2 - x1)
end

# left discrete Legendre transform for a free rigid body
function DLT2_linear(h, m, x1, x2)
	1.0 / h * m * (x2 - x1)
end

# discrete Newton for a free rigid body
function translational_integrator(h, m, x1, x2, x3, gravity, F) 
	DLT2_linear(h, m, x1, x2) + DLT1_linear(h, m, x2, x3) + h * m * gravity - F
end

# right discrete Legendre transform for a free rigid body
function DLT1_angular(h, J, q1, q2)
    H = [zeros(1, 3); I(3)]
    T = Diagonal([1.0; -1; -1; -1]) 
    (-4.0 / h) * transpose(L_mult(q1) * H) * T * transpose(R_mult(q2)) * H * J * transpose(H) * transpose(L_mult(q1)) * q2
end

# left discrete Legendre transform for a free rigid body
function DLT2_angular(h, J, q1, q2)
    H = [zeros(1, 3); I(3)]    
    (-4.0 / h) * transpose(L_mult(q2) * H) * L_mult(q1) * H * J * transpose(H) * transpose(L_mult(q1)) * q2
end

# discrete Euler-Lagrange for a free rigid body
function rotational_integrator(h, J, q1, q2, q3, τ)   
    DLT2_angular(h, J, q1, q2) + DLT1_angular(h, J, q2, q3) - 2.0 * τ
end

# finite-difference angular velocity 
function angular_velocity(h, q1, q2) 
    H = [zeros(1, 3); I(3)]
    2.0 * transpose(H) * transpose(L_mult(q1)) * (q2 - q1) / h
end

# Cayley map (modified-Rodriques parameter to quaternion)
function cayley(x)
    1.0 / sqrt(1.0 + norm(x)^2.0) * [1.0; x]
end

function quaternion_map(ω, h) 
	[sqrt(4 / h^2 - dot(ω, ω)); ω]
end

function φ(ϕ)
	1.0 / sqrt(1.0 + norm(ϕ)^2.0) * [1.0; ϕ]
end

kinematics(x, q, p) = x + rotation_matrix(q) * p 

function set_velocity(xa, va, qa, ωa, pa, xb, vb, qb, ωb, pb; Δv=[0.0; 0.0; 0.0], Δω=[0.0; 0.0; 0.0]) 

    # Ω(B/W)b = Ra->b * [Ω(B/A)a + Ω(A/W)a]
	ωb = rotation_matrix(L_mult(inverse_quaternion(qb)) * qa) * (Δω + ωa)
    ω1w = rotation_matrix(qa) * ωa
    ω2w = rotation_matrix(qb) * ωb
    Δvw = rotation_matrix(qa) * Δv
    cApB_w = (xb + rotation_matrix(qb) * pb) - xa
    pBcB_w = - rotation_matrix(qb) * pb
    vb = copy(va)
    vb += skew(ω1w) * cApB_w
    vb += skew(ω2w) * pBcB_w
    vb += Δvw

	return xb, vb, qb, ωb
end

function position_constraint(xa, qa, pa, xb, qb, pb) 
	ka = kinematics(xa, qa, pa) 
	kb = kinematics(xb, qb, pb)
	e = kb - ka 
	rotation_matrix(inverse_quaternion(qa)) * e 
end

using Symbolics 
@variables xa[1:3] qa[1:4] pa[1:3] xb[1:3] qb[1:4] pb[1:3] 

pc = position_constraint(xa, qa, pa, xb, qb, pb)

pc_xa = Symbolics.jacobian(pc, xa)
pc_qa = Symbolics.jacobian(pc, qa) * attitude_jacobian(qa)
pc_xb = Symbolics.jacobian(pc, xb)
pc_qb = Symbolics.jacobian(pc, qb) * attitude_jacobian(qb)

position_constraint_jacobian_parent_position = eval(Symbolics.build_function(pc_xa, xa, qa, pa, xb, qb, pb)[1])
position_constraint_jacobian_parent_orientation = eval(Symbolics.build_function(pc_qa, xa, qa, pa, xb, qb, pb)[1])
position_constraint_jacobian_child_position = eval(Symbolics.build_function(pc_xb, xa, qa, pa, xb, qb, pb)[1])
position_constraint_jacobian_child_orientation = eval(Symbolics.build_function(pc_qb, xa, qa, pa, xb, qb, pb)[1])

function visualize!(vis, p, q; timestep = 0.1)
	setvisible!(vis["/Background"], true)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)

    setobject!(vis[:pbody],
    	Rect(Vec(-0.5, -0.1, -0.1),Vec(1.0, 0.2, 0.2)),
    	MeshPhongMaterial(color = RGBA(1.0, 0.0, 0.0, 1.0)))

	setobject!(vis[:cbody],
		Rect(Vec(-0.5, -0.1, -0.1),Vec(1.0, 0.2, 0.2)),
		MeshPhongMaterial(color = RGBA(0.0, 1.0, 0.0, 1.0)))

   
    anim = MeshCat.Animation(convert(Int, floor(1.0 / timestep)))

     # for t = 1:length(q)
    	#  MeshCat.atframe(anim, t) do
    	# 	 settransform!(vis["satellite"],
    	# 		   compose(Translation((q[t][1:3] + [-0.25; -0.25; -0.25])...),
    	# 				 LinearMap(UnitQuaternion(q[t][4:7]...))))
    	#  end
     # end

	 for t = 1:length(q)
		xi = [q[t][(i-1) * 7 .+ (1:3)] for i = 1:N] 
		qi = [q[t][(i-1) * 7 + 3 .+ (1:4)] for i = 1:N] 
		MeshCat.atframe(anim, t) do
			for i = 1:N
				MeshCat.settransform!(vis[Symbol("body$i")],
					MeshCat.compose(MeshCat.Translation(xi[i]...),
							MeshCat.LinearMap(UnitQuaternion(qi[i]...))))
			end
		end
	 end
	#
	#
    MeshCat.setanimation!(vis, anim)
end

vis = Visualizer()
render(vis)

visualize!(vis, nothing, [z1, z2]; timestep=h)

# system parameters
gravity = [0.0; 0.0; 0.0 * -9.81]

ma = 1.0 
Ja = Diagonal([1.0; 1.0; 1.0])
pa = [0.5; 0.0; 0.0]

mb = 1.0 
Jb = Diagonal([1.0; 1.0; 1.0])
pb = [-0.5; 0.0; 0.0]

xa2 = [-0.5; 0.0; 0.0]
qa2 = [1.0; 0.0; 0.0; 0.0] 
va15 = [0.0; 0.0; 0.0]
ωa15 = [0.0; 0.0; 0.0] 
xa1 = xa2 - va15 * h 
qa1 = L_mult(qa2) * quaternion_map(-ωa15, h) * h / 2

xb2 = [0.5; 0.0; 0.0] 
qb2 = [1.0; 0.0; 0.0; 0.0]
vb15 = [0.0; 0.0; 0.0]
ωb15 = [0.0; 0.0; 0.0] 

xb2, vb15, qb2, ωb15 = set_minimal_velocity(xa2, va15, qa2, ωa15, pa, xb2, vb15, qb2, ωb15, pb, Δv=[0.0; 0.0; 0.0], Δω=[1.0; 0.0; 1.0])

xb1 = xb2 - vb15 * h 
qb1 = L_mult(qb2) * quaternion_map(-ωb15, h) * h / 2

ka1 = kinematics(xa1, qa1, pa)
kb1 = kinematics(xb1, qb1, pb) 

ka2 = kinematics(xa2, qa2, pa)
kb2 = kinematics(xb2, qb2, pb) 

mass = [ma, mb] 
inertia = [Ja, Jb] 
point = [pa, pb]

N = 2
M = 1
nz = N * (3 + 4)
ny = M * 3

function integrator(z1, z2, z3, y) 
	x1 = [z1[(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q1 = [z1[(i-1) * 7 + 3 .+ (1:4)] for i = 1:N] 
	x2 = [z2[(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q2 = [z2[(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]
	x3 = [z3[(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q3 = [z3[(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]

	[
		translational_integrator(h, ma, x1[1], x2[1], x3[1], gravity, zeros(3)) - transpose(position_constraint_jacobian_parent_position(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y;
		translational_integrator(h, mb, x1[2], x2[2], x3[2], gravity, zeros(3)) - transpose(0.5 * position_constraint_jacobian_child_position(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y;
		rotational_integrator(h, Ja, q1[1], q2[1], q3[1], zeros(3)) - transpose(position_constraint_jacobian_parent_orientation(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y;
		rotational_integrator(h, Jb, q1[2], q2[2], q3[2], zeros(3)) - transpose(0.5 * position_constraint_jacobian_child_orientation(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y;
		position_constraint(x3[1], q3[1], pa, x3[2], q3[2], pb);
	]
end

function system_attitude_jacobian(z) 
	q = [z[(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]
	cat([cat(I(3), attitude_jacobian(q[i]), dims=(1,2)) for i = 1:N]..., I(3), dims=(1,2))
end

z1 = [xa1; qa1; xb1; qb1]
z2 = [xa2; qa2; xb2; qb2]
z3 = copy(z2)
y = zeros(ny)

# integrator(z1, z2, z3, y) 
# system_attitude_jacobian(z3)

h = 0.01
T = 250

# variational
function newton_variational(z1, z2)
	y = zeros(ny) 
	z = copy(z2) 
	w = [z; y]

	_r(a) = integrator(z1, z2, a[1:nz], a[nz .+ (1:ny)])

	function cand(x, Δ, α)
		x⁺ = copy(x)
		# body a
		x⁺[1:3] = x[1:3] + α * Δ[1:3] 
		x⁺[3 .+ (1:4)] = L_mult(x[3 .+ (1:4)]) * φ(α * Δ[3 .+ (1:3)])

		# body b 
		x⁺[7 .+ (1:3)] = x[7 .+ (1:3)] + α * Δ[6 .+ (1:3)] 
		x⁺[7 + 3 .+ (1:4)] = L_mult(x[7 + 3 .+ (1:4)]) * φ(α * Δ[6 + 3 .+ (1:3)])

		# impulse 
		x⁺[7 + 3 + 4 .+ (1:3)] = x[7 + 3 + 4 .+ (1:3)] + α * Δ[6 + 3 + 3 .+ (1:3)]

		return x⁺
	end

	r = _r(w)
	# @show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, w) * system_attitude_jacobian(w)
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		w_cand = cand(w, Δ, α)
		r_cand = _r(w_cand)
		# @show norm(r_cand)
		iter = 0

		while norm(r_cand) > norm(r)
			α *= 0.5
			w_cand = cand(w, Δ, α)
			r_cand = _r(w_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		w = w_cand
		r = r_cand

		if norm(r) < 1.0e-10
			return w
		end
	end
	@warn "newton failure"
	return w
end


# newton_variational(z1, z2)

z = [z1, z2]

H_linear_a = []#DLT2_linear(h, ma, z[end-1][1:3], q[end][1:3])]
H_angular_a = []#rotation_matrix(qa2) * DLT2_angular(h, Ja, z[end-1][3 .+ (1:4)], z[end][3 .+ (1:4)])]
H_linear_b = []#DLT2_linear(h, ma, z[end-1][7 .+ (1:3)], z[end][7 .+ (1:3)])]
H_angular_b = []#rotation_matrix(qb2) * DLT2_angular(h, Jb, z[end-1][10 .+ (1:4)], z[end][10 .+ (1:4)])]

H_linear = [] 
H_angular = [] 

for t = 1:T
	sol = newton_variational(z[end-1], z[end])
	push!(z, sol[1:nz])

	## momentum 
	x1 = [z[end-2][(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q1 = [z[end-2][(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]

	x2 = [z[end-1][(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q2 = [z[end-1][(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]

	x3 = [z[end][(i-1) * 7 .+ (1:3)] for i = 1:N] 
	q3 = [z[end][(i-1) * 7 + 3 .+ (1:4)] for i = 1:N]

	y = sol[nz .+ (1:ny)] 

	push!(H_linear_a, DLT1_linear(h, ma, x2[1], x3[1]) - 0.5 * transpose(position_constraint_jacobian_parent_position(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y)
	push!(H_angular_a, rotation_matrix(q2[1]) * (DLT1_angular(h, Ja, q2[1], q3[1]) - 0.25 * transpose(position_constraint_jacobian_parent_orientation(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y))
	push!(H_linear_b, DLT1_linear(h, mb, x2[2], x3[2]) - 0.5 * transpose(position_constraint_jacobian_child_position(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y)
	push!(H_angular_b, rotation_matrix(q2[2]) * (DLT1_angular(h, Jb, q2[2], q3[2]) - 0.25 * transpose(position_constraint_jacobian_child_orientation(x2[1], q2[1], pa, x2[2], q2[2], pb)) * y))

	@show com = (ma * x2[1] + mb * x2[2]) / (ma + mb) 

	linear_total = H_linear_a[end] + H_linear_b[end]
	v_com = linear_total ./ (ma + mb) 

	angular_total = H_angular_a[end] + H_angular_b[end] 
	angular_total += cross(x2[1] - com, ma .* (H_linear_a[end] ./ ma - linear_total ./ (ma + mb)))
	angular_total += cross(x2[2] - com, mb .* (H_linear_b[end] ./ mb - linear_total ./ (ma + mb)))
	
	push!(H_linear, linear_total) 
	push!(H_angular, angular_total)
end

H_linear
H_angular

visualize!(vis, nothing, z, timestep = h)

# plot(hcat(H...)')
