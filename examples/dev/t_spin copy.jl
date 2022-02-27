# skew-symmetric matrix
function hat(x)
    return [0 -x[3] x[2];
            x[3] 0 -x[1];
           -x[2] x[1] 0]
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
function DLT1(h, J, q1, q2)
    H = [zeros(1, 3); I(3)]
    T = Diagonal([1.0; -1; -1; -1]) 
    (-4.0 / h) * transpose(L_mult(q1) * H) * T * transpose(R_mult(q2)) * H * J * transpose(H) * transpose(L_mult(q1)) * q2
end

# left discrete Legendre transform for a free rigid body
function DLT2(h, J, q1, q2)
    H = [zeros(1, 3); I(3)]    
    (-4.0 / h) * transpose(L_mult(q2) * H) * L_mult(q1) * H * J * transpose(H) * transpose(L_mult(q1)) * q2
end

# discrete Euler-Lagrange equation for a free rigid body
function integrator(h, J, q1, q2, q3, τ)   
    DLT2(h, J, q1, q2) + DLT1(h, J, q2, q3) + 2.0 * τ
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

function visualize!(vis, p, q; timestep=0.1)
	setvisible!(vis["/Background"], true)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)

    setobject!(vis[:satellite],
    	Rect(Vec(-0.25, -0.25, -0.25),Vec(0.5, 0.5, 0.5)),
    	MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    arrow_x = ArrowVisualizer(vis[:satellite][:arrow_x])
    mat = MeshPhongMaterial(color=RGBA(1.0, 0.0, 0.0, 1.0))
    setobject!(arrow_x, mat)
    settransform!(arrow_x,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.75, 0.0, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_y = ArrowVisualizer(vis[:satellite][:arrow_y])
    mat = MeshPhongMaterial(color=RGBA(0.0, 1.0, 0.0, 1.0))
    setobject!(arrow_y, mat)
    settransform!(arrow_y,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.75, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_z = ArrowVisualizer(vis[:satellite][:arrow_z])
    mat = MeshPhongMaterial(color=RGBA(0.0, 0.0, 1.0, 1.0))
    setobject!(arrow_z, mat)
    settransform!(arrow_z,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.0, 0.75),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    anim = MeshCat.Animation(convert(Int, floor(1.0 / timestep)))

     # for t = 1:length(q)
    	#  MeshCat.atframe(anim, t) do
    	# 	 settransform!(vis["satellite"],
    	# 		   compose(Translation((q[t][1:3] + [-0.25; -0.25; -0.25])...),
    	# 				 LinearMap(UnitQuaternion(q[t][4:7]...))))
    	#  end
     # end

	 for t = 1:length(q)
		MeshCat.atframe(anim, t) do
			MeshCat.settransform!(vis["satellite"],
				  MeshCat.compose(MeshCat.Translation(([-0.0; -0.0; -0.0])...),
						MeshCat.LinearMap(UnitQuaternion(q[t][1:4]...))))
		end
	 end
	#
	#
    MeshCat.setanimation!(vis, anim)
end

vis=visualizer()
render(vis)

nq = 4
# nv = 3
H = [zeros(1, 3); I]

J = Diagonal([1.0; 2.0; 3.0])

h = 0.01
T = 1000

# dynamics
function f(x, u)
	q = x[1:4]
	ω = x[5:7]
	τ = u[1:3]

	[0.5 * L_mult(q) * H * ω;
	 J \ (τ - cross(ω, J * ω))]
end

# explicit
function f(x, u, h)
    x + h * f(x + 0.5 * h * f(x, u), u) # explicit midpoint
end

# implicit
function f(x⁺, x, u, h)
	x⁺ - (x + h * f(0.5 * (x + x⁺), u)) # implicit midpoint
end


function angular_momentum(x) 
    q = x[1:4] 
    ω = x[5:7] 
    rotation_matrix(q) * J * ω 
end

q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 10.0; 0.1]

x1 = [q1; ω1]
u1 = zeros(3)

function newton(x, u, h)
	x⁺ = copy(x)
	_r(z) = f(z, x, u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		x⁺[1:4] = L_mult(x[1:4]) * φ(α * Δ[1:3])
		x⁺[5:7] += α * Δ[4:6]

		return x⁺
	end

	r = _r(x⁺)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, x⁺) * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		x_cand = cand(x⁺, Δ, α)
		r_cand = _r(x_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(x⁺, Δ, α)
			r_cand = _r(x_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		x⁺ = x_cand
		r = r_cand

		if norm(r) < 1.0e-8
			return x⁺
		end
	end
	@warn "newton failure"
	return x
end
newton(x1, u1, h)

x = [copy(x1)]
for t = 1:T
	push!(x, newton(x[end], u1, h))
end

visualize!(vis, nothing, x, timestep= h)

am = [angular_momentum(xt) - angular_momentum(x[1]) for xt in x]
plot(hcat(am...)')

# variational
function newton_variational(q1, q2, u, h)
	q3 = copy(q2) 

	_r(z) = integrator(h, J, q1, q2, z, u)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		x⁺[1:4] = L_mult(x[1:4]) * φ(α * Δ[1:3])
		# x⁺[1:3] += α * Δ[1:3]
		return x⁺
	end

	r = _r(q3)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, q3) * attitude_jacobian(q3)
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		q_cand = cand(q3, Δ, α)
		r_cand = _r(q_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) > norm(r)
			α *= 0.5
			q_cand = cand(q3, Δ, α)
			r_cand = _r(q_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		q3 = q_cand
		r = r_cand

		if norm(r) < 1.0e-10
			return q3
		end
	end
	@warn "newton failure"
	return q3
end

q2 = [1.0; 0.0; 0.0; 0.0]
ω15 = [0.0; 10.0; 0.01] 
q1 = L_mult(q2) * quaternion_map(-ω15, h) * h / 2
u1 = zeros(3) 

newton_variational(q1, q2, u1, h)

q = [q1, q2]
H = [rotation_matrix(q[end]) * DLT2(h, J, q[end-1], q[end])]
for t = 1:T
	push!(q, newton_variational(q[end-1], q[end], u1, h))
	push!(H, rotation_matrix(q[end]) * DLT2(h, J, q[end-1], q[end]))
end

visualize!(vis, nothing, q, timestep= h)

plot(hcat(H...)')
