function lineSearch!(mechanism::Mechanism, α, rvio, bvio, opts; warning::Bool = false)
    scale = 0
    system = mechanism.system
    eqcs = mechanism.eqconstraints
    bodies = mechanism.bodies
    ineqcs = mechanism.ineqconstraints

    rvio_cand, bvio_cand = Inf * ones(2)
    for n = Base.OneTo(opts.max_ls)
        for ineqc in ineqcs
            lineStep!(α, ineqc, getentry(system, ineqc.id), scale)
        end
        for eqc in eqcs
            lineStep!(α, eqc, getentry(system, eqc.id), scale)
        end
        for body in bodies
            ϕmax = 3.9 / mechanism.Δt^2
            lineStep!(α, body, getentry(system, body.id), scale, ϕmax = ϕmax)
            if dot(body.state.ϕsol[2], body.state.ϕsol[2]) > 3.91 / mechanism.Δt^2
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

    warning && (@info string("lineSearch! did not converge. n = ", iter, ". Last tol: ", meritf1))
    return rvio_cand, bvio_cand
end

@inline function lineStep!(α::T, body::Body, vector_entry::Entry, scale; ϕmax = Inf) where T
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[1; 2; 3]]
    body.state.ϕsol[2] = body.state.ϕsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[4; 5; 6]]
    # @warn "clipping ϕ might not be the best solution, we might add an IP constraints on this variable instead."
    ϕ = body.state.ϕsol[2]
    ϕnorm = dot(ϕ, ϕ)
    if ϕnorm > ϕmax
        body.state.ϕsol[2] *= ϕmax / ϕnorm
    end
    return
end

@inline function lineStep!(α::T, eqc::EqualityConstraint, vector_entry::Entry, scale) where T
    eqc.λsol[2] = eqc.λsol[1] + 1.0 / (2^scale) * vector_entry.value
    return
end

@inline function lineStep!(α::T, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, scale) where {T,N,Nc,Cs,N½}
    ineqc.ssol[2] = ineqc.ssol[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(1:N½)]
    ineqc.γsol[2] = ineqc.γsol[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(N½+1:N)]
    return
end
