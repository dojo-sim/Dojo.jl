function line_search!(mechanism::Mechanism, α, rvio, bvio, opts; warning::Bool = false)
    scale = 0
    system = mechanism.system

    rvio_cand, bvio_cand = Inf * ones(2)
    for n = Base.OneTo(opts.max_ls)
        for contact in mechanism.contacts
            candidate_step!(α, contact, get_entry(system, contact.id), scale)
        end
        for joint in mechanism.joints
            candidate_step!(α, joint, get_entry(system, joint.id), scale)
        end
        for body in mechanism.bodies
            ϕmax = 3.9 / mechanism.timestep^2
            candidate_step!(α, mechanism, body, get_entry(system, body.id), scale, ϕmax = ϕmax)
            if dot(body.state.ϕsol[2], body.state.ϕsol[2]) > 3.91 / mechanism.timestep^2
                error("Excessive angular velocity. Body-ID: $(string(body.name)) " * string(body.id) * ", ω: " * string(body.state.ϕsol[2]) * ".")
            end
        end

        rvio_cand = residual_violation(mechanism)
        bvio_cand = bilinear_violation(mechanism)

        if (rvio_cand > rvio) && (bvio_cand > bvio)
            scale += 1
        else
            return rvio_cand, bvio_cand
        end
    end

    warning && (@info string("line_search! did not converge. n = ", iter, ". Last tol: ", meritf1))
    return rvio_cand, bvio_cand
end

@inline function candidate_step!(α, mechanism::Mechanism, body::Body, vector_entry::Entry, scale; ϕmax = Inf)
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[1; 2; 3]]
    body.state.ϕsol[2] = body.state.ϕsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[4; 5; 6]]
    ϕ = body.state.ϕsol[2]
    ϕdot = dot(ϕ, ϕ)
    if ϕdot > ϕmax
        println("clipping ", scale, scn((ϕdot - ϕmax) / ϕmax), " ", scn(ϕdot), " ", scn(ϕmax), " ", body.name)
        body.state.ϕsol[2] *= ϕmax / ϕdot # this is overkill, but works better than sqrt(ϕmax/ϕdot)
    end
    return
end

@inline function candidate_step!(α, joint::JointConstraint, vector_entry::Entry, scale)
    joint.impulses[2] = joint.impulses[1] + 1.0 / (2^scale) * α * vector_entry.value
    return
end

@inline function candidate_step!(α::T, contact::ContactConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, scale) where {T,N,Nc,Cs,N½}
    contact.primal[2] = contact.primal[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(1:N½)]
    contact.dual[2] = contact.dual[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(N½+1:N)]
    return
end
