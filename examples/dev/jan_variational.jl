using LinearAlgebra, Plots, ForwardDiff
function qmult(q1, q2)
    v1 = q1[2:4]
    s1 = q1[1]
    
    v2 = q2[2:4]
    s2 = q2[1]
    
    [s1*s2 - v1'*v2; s1*v2 + s2*v1 + hat(v1)*v2]
end

function qrot(q, r)
    
    v = q[2:4]
    vh = hat(v);
    w = q[1];
    
    rrot = r + 2*vh*(vh*r + w*r); 
end

function hat(x)
    [  0   -x[3]  x[2]
            x[3]   0   -x[1]
        -x[2]  x[1]  0];
end

function skewplusdiag(v, w) where T
    [
         w    -v[3]  v[2]
         v[3]  w    -v[1]
        -v[2]  v[1]  w
    ]
end

function L_multiply(q)
	s = q[1]
	v = q[2:4]


	[s -transpose(v);
     v s * I + hat(v)]
end



function newton(x, p)
	x⁺ = copy(x)
	_r(phi) = p - (sqrt(4 / dt^2.0 - phi' * phi) * J * phi - hat(phi) * (J*phi))

	function cand(x, Δ, α)
		x⁺ = copy(x)
		x⁺[1:3] += α * Δ[1:3]
		return x⁺
	end

	r = _r(x⁺)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, x⁺)
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

		if norm(r) < 1.0e-12
			return x⁺
		end
	end
	@warn "newton failure"
	return x
end

q = [[1.0; 0.0; 0.0; 0.0],]
J = Diagonal([1.0; 2.0; 3.0])
w0 = [pi/10; pi/6; pi/8]
p = J*w0
w = [J * w0]
dt = 0.01
phi = w0#*dt/2
N = 1000

phi = newton(phi, p)
# phi = w0
for k = 1:N
    # gg = .5*dt*p
    p = (sqrt(4 / dt^2.0 - phi' * phi) * J * phi - hat(phi) * (J*phi))
    gg = p
    for j = 1:10
        sq2 = sqrt(4 / dt^2 - phi' * phi)
        e = (J * phi) * sqrt(4 / dt^2 - phi' * phi) + cross(phi, J * phi) - gg;
        dedphi = skewplusdiag(phi, sq2) * J - J * phi * (phi' / sq2) - hat(J * phi)
        phi = phi - dedphi \ e;
    end
    push!(q, 0.5 * dt * L_multiply(q[end]) * [sqrt(4.0 / dt^2.0 - phi' * phi); phi])

    p = (sqrt(4 / dt^2.0 - phi' * phi) * J * phi - hat(phi) * (J*phi))
    push!(w, p)
end

h = [qrot(q[t], w[t]) - qrot(q[1], w[1]) for t = 1:length(q)-1]
plot(hcat(h...)')