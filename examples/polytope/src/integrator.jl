################################################################################
# Methods
################################################################################
function implicit_integrator(fixed_collider, soft, timestep, z)
    opts = soft.options
	x2, v15, q2, ϕ15 = unpack_maximal_state(z, 1)
	x1 = Dojo.next_position(x2, -v15, timestep)
	q1 = Dojo.next_orientation(q2, -ϕ15, timestep)

    soft.x = x2
    soft.q = q2
    ψ, impact_normal, barycenter = collision(fixed_collider, soft)
    coulomb_direction(v) = - atan(opts.coulomb_smoothing * norm(v)) * v/(opts.coulomb_regularizer + norm(v))

    function residual(vϕ)
        v25 = vϕ[SUnitRange(1,3)]
        ϕ25 = vϕ[SUnitRange(4,6)]
        x3 = Dojo.next_position(x2, v25, timestep)
        q3 = Dojo.next_orientation(q2, ϕ25, timestep)

        # dynamics
        D1x = - 1.0 / timestep * mass * (x2 - x1) - 0.5 * timestep * mass * gravity
        D2x =   1.0 / timestep * mass * (x3 - x2) - 0.5 * timestep * mass * gravity
        D1q = -2.0 / timestep * LVᵀmat(q2)' * Lmat(q1) * Vᵀmat() * inertia * Vmat() * Lmat(q1)' * vector(q2)
        D2q = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)
        dynT = D2x + D1x
        dynR = D2q + D1q
        d = [dynT; dynR]

        # inputs
        barycenter_w = Dojo.vector_rotate(barycenter, q2) # TODO could be q3
        vc = v25 + Dojo.vector_rotate(Dojo.skew(-barycenter) * ϕ25, q2) # TODO could be q3
        vc_normal = vc' * impact_normal * impact_normal
        vc_tangential = vc - vc_normal

        F_impact = -opts.impact_damper * ψ * vc_normal
        F_impact += opts.impact_spring * ψ * impact_normal
		F_friction = -opts.sliding_drag * norm(F_impact) * vc_tangential
		F_friction += opts.sliding_friction * norm(F_impact) * coulomb_direction(vc_tangential)
        F_contact = F_impact + F_friction
        F2 = mass*gravity + F_contact

        τ_w = -Dojo.skew(F_contact) * barycenter_w
        τ2 = Dojo.vector_rotate(τ_w, inv(q2)) # TODO could be q3
        τ2 += -opts.rolling_drag * norm(F_impact) * ϕ25
        τ2 += opts.rolling_friction * norm(F_impact) * coulomb_direction(ϕ25)

		d -= timestep * [F2; τ2]
        return d
    end

    vϕ = newton_solver(residual)
	v25 = vϕ[SUnitRange(1,3)]
	ϕ25 = vϕ[SUnitRange(4,6)]
	x3 = Dojo.next_position(x2, v25, timestep)
	q3 = Dojo.next_orientation(q2, ϕ25, timestep)
    z1 = [x3; v25; Dojo.vector(q3); ϕ25]
    return z1, ψ
end

function newton_solver(residual)
	x = zeros(6)
	for i = 1:10
		res = residual(x)
		(norm(res, Inf) < 1e-6) && break
		jac = FiniteDiff.finite_difference_jacobian(x -> residual(x), x)
		Δ = - jac \ res
		α = inegrator_line_search(x, res, Δ, residual)
		x += α * Δ
	end
	return x
end

function inegrator_line_search(x, res, Δ, residual)
    α = 1.0
    for i = 1:10
        (norm(residual(x + α * Δ), Inf) <= norm(res, Inf)) && break
        α *= 0.5
    end
    return α
end

function implicit_simulation(fixed_collider, soft, timestep, N, z0)
    z = [z0]
    ψ = []
    for i = 1:N
        @show i
        z0, ψ0 = implicit_integrator(fixed_collider, soft, timestep, z0)
        push!(z, z0)
        push!(ψ, ψ0)
    end
    return z, ψ
end
