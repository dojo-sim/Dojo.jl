function conjugate(q)
	s = q[1]
	v = q[2:4]

	return [s; -v]
end

function L_multiply(q)
	s = q[1]
	v = q[2:4]


	SMatrix{4,4}([s -transpose(v);
	              v s * I + skew(v)])
end

function R_multiply(q)
	s = q[1]
	v = q[2:4]


	SMatrix{4,4}([s -transpose(v);
	              v s * I - skew(v)])
end

function multiply(q1, q2)
	L_multiply(q1) * q2
end

function qmult(q1, q2)
    v1 = q1[1:3]
    s1 = q1[4]
    
    v2 = q2[1:3]
    s2 = q2[4]
    
    [s1*v2 + s2*v1 + skew(v1)*v2; s1*s2 - v1'*v2]
end

# eq. 14 http://roboticexplorationlab.org/papers/planning_with_attitude.pdf
function attitude_jacobian(q)
	s = q[1]
	v = q[2:4]

	[-transpose(v);
	 s * I + skew(v)]
end

# eq. 16 http://roboticexplorationlab.org/papers/maximal_coordinate_dynamics.pdf
function ω_finite_difference(q1, q2, h)
	2.0 * multiply(conjugate(q1), (q2 - q1) ./ h)[2:4]
end

# https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
function quaternion_rotation_matrix(q)
	r, i, j, k  = q

	r11 = 1.0 - 2.0 * (j^2.0 + k^2.0)
	r12 = 2.0 * (i * j - k * r)
	r13 = 2.0 * (i * k + j * r)

	r21 = 2.0 * (i * j + k * r)
	r22 = 1.0 - 2.0 * (i^2.0 + k^2.0)
	r23 = 2.0 * (j * k - i * r)

	r31 = 2.0 * (i * k - j * r)
	r32 = 2.0 * (j * k + i * r)
	r33 = 1.0 - 2.0 * (i^2.0 + j^2.0)

	SMatrix{3,3}([r11 r12 r13;
	              r21 r22 r23;
				  r31 r32 r33])
end

# Cayley map
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

function angular_momentum(x) 
    q = x[1:4] 
    ω = x[5:7] 
    quaternion_rotation_matrix(q) * J * ω 
end

vis=visualizer()
render(vis)

nq = 4
# nv = 3
H = [zeros(1, 3); I]

J = Diagonal([1.0; 1.0; 1.0])

h = 0.01
T = 1000

function f(x, u)
	q = x[1:4]
	ω = x[5:7]
	τ = u[1:3]

	[0.5 * L_multiply(q) * H * ω;
	 J \ (τ - cross(ω, J * ω))]
 end

function f(x, u, h)
    x + h * f(x + 0.5 * h * f(x, u), u) # explicit midpoint
end

q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 10.0; 0.1]

x1 = [q1; ω1]
u1 = zeros(3)

x = [copy(x1)]
for t = 1:T
	push!(x, f(x[end], u1, h))
end

visualize!(vis, nothing, x, timestep= h)

am = [angular_momentum(xt) - angular_momentum(x[1]) for xt in x]
plot(hcat(am...)')

# implicit
function f(x⁺, x, u, h)
	x⁺ - (x + h * f(0.5 * (x + x⁺), u)) # implciit midpoint
end

function newton(x, u, h)
	x⁺ = copy(x)
	_r(z) = f(z, x, u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
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

function newton_variational(x, u, h)
	q2 = copy(x)[1:4]
	ω2 = copy(x)[5:7]

	function f(ω1, ω2, q2, u, h)
		(J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
			+ cross(ω1, J * ω1)
			- 2.0 * u[1:3])
	end

	_r(z) = f(x[5:7], z, x[1:4], u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		# x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
		x⁺[1:3] += α * Δ[1:3]

		return x⁺
	end

	r = _r(ω2)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, ω2)# * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		ω_cand = cand(ω2, Δ, α)
		r_cand = _r(ω_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(ω2, Δ, α)
			r_cand = _r(ω_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		ω2 = ω_cand
		r = r_cand

		if norm(r) < 1.0e-10
			q3 = 0.5 * h * L_multiply(q2) * [sqrt((2.0 / h)^2.0 - ω2' * ω2); ω2]
			return [q3; ω2]
		end
	end
	@warn "newton failure"
	return x
end

newton_variational(x1, u1, h)

x = [copy(x1)]
for t = 1:T
	push!(x, newton_variational(x[end], u1, h))
end

visualize!(vis, nothing, x, timestep= h)

function newton_zac_variational(x, u, h; p = zeros(3), init=false)
	q2 = copy(x)[1:4]
	ω2 = copy(x)[5:7]

	function f(ω1, ω2, q2, u, h)
        if init 
		    sqrt(1.0 - ω2' * ω2) * J * ω2 + cross(ω2, J * ω2) - 0.5 * h * p + 0.5 * h^2.0 * u
        else
		    sqrt(1.0 - ω2' * ω2) * J * ω2 + cross(ω2, J * ω2) - (sqrt(1.0 - ω1' * ω1) * J * ω1 - cross(ω1, J * ω1)) + 0.5 * h^2.0 * u
        end
	end

	_r(z) = f(x[5:7], z, x[1:4], u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		# x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
		x⁺[1:3] += α * Δ[1:3]

		return x⁺
	end

	r = _r(ω2)
	@show norm(r)

	for i = 1:25
		∇r = ForwardDiff.jacobian(_r, ω2)# * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		ω_cand = cand(ω2, Δ, α)
		r_cand = _r(ω_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(ω2, Δ, α)
			r_cand = _r(ω_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		ω2 = ω_cand
		r = r_cand

		if norm(r) < 1.0e-8
			q3 = qmult(q2, [ω2; sqrt(1.0 - ω2' * ω2)])
			return [q3; ω2]
		end
	end
	@warn "newton failure"
	return x
end

function newton_zac_variational_init(x, u, h; p = zeros(3), init=false)
	q2 = copy(x)[1:4]
	ω2 = copy(x)[5:7]

	function f(ω1, ω2, q2, u, h)
        if init 
		    sqrt(1.0 - ω2' * ω2) * J * ω2 + cross(ω2, J * ω2) - 0.5 * h * p + 0.5 * h^2.0 * u
        else
		    sqrt(1.0 - ω2' * ω2) * J * ω2 + cross(ω2, J * ω2) - (sqrt(1.0 - ω1' * ω1) * J * ω1 - cross(ω1, J * ω1)) + 0.5 * h^2.0 * u
        end
	end

	_r(z) = f(x[5:7], z, x[1:4], u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		# x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
		x⁺[1:3] += α * Δ[1:3]

		return x⁺
	end

	r = _r(ω2)
	@show norm(r)

	for i = 1:25
		∇r = ForwardDiff.jacobian(_r, ω2)# * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		ω_cand = cand(ω2, Δ, α)
		r_cand = _r(ω_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(ω2, Δ, α)
			r_cand = _r(ω_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		ω2 = ω_cand
		r = r_cand

		if norm(r) < 1.0e-8
			q3 = qmult(q2, [ω2; sqrt(1.0 - ω2' * ω2)])
			return [q3; ω2]
		end
	end
	@warn "newton failure"
	return x
end

q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 1.0; 0.01]
x1 = [q1; 0.5 * h * ω1]
p1 = J * ω1
newton_zac_variational(x1, u1, h)
newton_zac_variational(x1, u1, h, contact_point=contact_point1, init=true)

x = [copy(x1)]
for t = 1:T
    if t == 1
        push!(x, newton_zac_variational(x[end], u1, h, contact_point=contact_point1, init=true))
    else
	    push!(x, newton_zac_variational(x[end], u1, h))
    end
end

vis=visualizer() 
render(vis)
visualize!(vis, nothing, x, timestep= h)

function angular_momentum_variational(x) 
    q = x[1:4] 
    ϕ = x[5:7] 
    # quaternion_rotation_matrix(q) * (2.0 / h * sqrt(1.0 - ϕ' * ϕ) * J * ϕ - cross(ϕ, J * ϕ))
    quaternion_rotation_matrix(q) * (2.0 / h * sqrt(1.0 - ϕ' * ϕ) * J * ϕ + cross(ϕ, J * ϕ))
end

am = [Array(angular_momentum_variational(xt) - angular_momentum_variational(x[1])) for xt in x]
plot(hcat(am...)')

0.5 * h * p1

function newton_variational_2(q1, q2, u, h)
	q3 = copy(q2)

	# @show ω_finite_difference(q1, q2, h)
	# @show norm(q3)
	# @show q3
	function f(q1, q2, q3, u, h)
		ω1 = ω_finite_difference(q1, q2, h)
		ω2 = ω_finite_difference(q2, q3, h)

		(J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
			+ cross(ω1, J * ω1)
			- 2.0 * u[1:3])
	end

	_r(z) = f(q1, q2, z, u, h)

	function cand(q3, Δ, α)
		return L_multiply(q3) * φ(α * Δ[1:3])
	end

	r = _r(q3)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, q3) * attitude_jacobian(q3)
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		q3_cand = cand(q3, Δ, α)
		r_cand = _r(q3_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) > norm(r)
			α *= 0.5
			q3_cand = cand(q3, Δ, α)
			r_cand = _r(q3_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		q3 = q3_cand
		r = r_cand

		# @show ω_finite_difference(q2, q3, h)
		# @show norm(q3)
		# @show q3
		if norm(r) < 1.0e-10
			return q3
		end
	end
	@warn "newton failure"
	return q2
end

x = [copy(q1), copy(Array(0.5 * h * L_multiply(q1) * [sqrt((2.0 / h)^2.0 - ω1' * ω1); ω1]))]

for t = 1:500
	push!(x, newton_variational_2(x[end-1], x[end], u, h))
end

visualize!(vis, nothing, x, timestep= h)
