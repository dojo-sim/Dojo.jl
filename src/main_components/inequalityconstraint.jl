
mutable struct InequalityConstraint{T,N,Nc,Cs,N½} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    # Currently only single constraint and child
    constraints::Cs # can be of type
    parentid::Int64
    childids::SVector{1,Union{Int64,Nothing}}
    ssol::Vector{SVector{N½,T}} # holds the slack variable
    γsol::Vector{SVector{N½,T}} # holds the dual of the slack variable

    function InequalityConstraint(data; name::String="")
        bound, parentid, childid = data
        T = getT(bound)

        childids = [childid]
        constraint = Tuple([bound])
        N = length(constraint[1])
        N½ = Int64(N/2)

        ssol = [neutral_vector(bound) for i = 1:2]
        γsol = [neutral_vector(bound) for i = 1:2]
        new{T,N,1,typeof(constraint),N½}(getGlobalID(), name, constraint, parentid, childids, ssol, γsol)
    end
end


function resetVars!(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    ineqc.ssol[1] = neutral_vector(ineqc.constraints[1])
    ineqc.ssol[2] = neutral_vector(ineqc.constraints[1])
    ineqc.γsol[1] = neutral_vector(ineqc.constraints[1])
    ineqc.γsol[2] = neutral_vector(ineqc.constraints[1])
    return
end

# contribution of the inequality constraint (impact or friction) to the dynamics equation d = 0
@inline function constraintForceMapping!(mechanism, body::Body, ineqc::InequalityConstraint)
    body.state.d -= ∂g∂ʳpos(mechanism, ineqc, body)' * ineqc.γsol[2]
    return
end

# contribution of the inequality constraint (impact or friction) to the dynamics equation d = 0
@inline function ∂constraintForceMapping!(mechanism, body::Body, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Δt = mechanism.Δt
    x3, q3 = posargsnext(body.state, Δt)
    x2, v2, q2, ω2 = fullargssol(body.state)

    for i=1:Nc
        bnd = ineqc.constraints[i]
        bnd_type = typeof(ineqc.constraints[i])

        M = ∂integration(q2, ω2, Δt)
        function d(vars)
            x = vars[1:3]
            q = UnitQuaternion(vars[4:7]..., false)
            return ∂g∂ʳpos(bnd, x, q)' * ineqc.γsol[2]
        end

        if bnd_type <: ContactBound
            body.state.D -= _dN(x3, [q3.w; q3.x; q3.y; q3.z], ineqc.γsol[2][1:1], bnd.p) * M
            body.state.D -= _dB(x3, [q3.w; q3.x; q3.y; q3.z], ineqc.γsol[2][2:4], bnd.p) * M
        elseif bnd_type <: ImpactBound11
            body.state.D -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
        elseif bnd_type <: LinearContactBound11
            body.state.D -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
        end
    end
    return
end


@inline function ∂gab∂ʳba(mechanism, body::Body, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    # derivative of gs and complementarity constraints wrt pos and vel
    Z = szeros(T,N½,6)
    # derivative of complementarity constraints wrt to pos and vel are always zero because it's independent of pos and vel.
    return [Z;-∂g∂ʳpos(mechanism, ineqc, body)]', [Z;∂g∂ʳvel(mechanism, ineqc, body)]
end
@inline function ∂gab∂ʳba(mechanism, ineqc1::InequalityConstraint, ineqc2::InequalityConstraint)
    G1, G2 = ∂gab∂ʳba(ineqc1.constraints[1], ineqc2.constraints[1])
    return G1, G2
end

function ∂g∂ʳposa(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳposa(ineqc.constraints[1], body, nothing, mechanism.Δt)
end

function ∂g∂ʳvela(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳvela(ineqc.constraints[1], body, nothing, mechanism.Δt)
end
