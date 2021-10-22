
# TODO currently just single parent and child
mutable struct InequalityConstraint{T,N,Nc,Cs,N½} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    # Currently only single constraint and child
    constraints::Cs # can be of type
        # ImpactBound
        # FrictionBound
        # BetaBound
        # ConeBound
    parentid::Int64
        # parent = body for impact bound
        # parent = friction for friction bound
        # parent = friction for beta bound
    childids::SVector{1,Union{Int64,Nothing}}
        # nothing for impact bound
        # impact bound for friction bound
        # nothing for beta bound
    ssol::Vector{SVector{N½,T}} # holds the slack variable
        # sψ for FrictionBound
        # sη for BetaBound
        # sγ for ImpactBound
    γsol::Vector{SVector{N½,T}} # holds the dual of the slack variable
        # ψ for FrictionBound
        # η for BetaBound
        # γ for ImpactBound

    function InequalityConstraint(data; name::String="")
        bound, parentid, childid = data
        T = getT(bound)

        childids = [childid]
        constraint = Tuple([bound])
        N = length(constraint[1])
        # @show N
        # @show typeof(constraint[1])
        N½ = Int64(N/2)

        # TODO why is there two elements in his vector? (linesearch probably)
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
        if typeof(ineqc.constraints[i]) <: ContactBound
            # @show "bound!"
            M = [Δt * I zeros(3,3); zeros(4,3) Lmat(q2)*derivωbar(ω2, Δt)*Δt/2]
            body.state.D -= _dN(x3, [q3.w; q3.x; q3.y; q3.z], ineqc.γsol[2][1:1], bnd.p) * M
            body.state.D -= _dB(x3, [q3.w; q3.x; q3.y; q3.z], ineqc.γsol[2][2:4], bnd.p) * M
        end
    end
    return
end

function g(mechanism, ineqc::InequalityConstraint)
    # this is the incomplete residual (we are missing -1 * slacks) of the equality constraint in the friction problem
    # specialized function in friction and beta_bounds
    g(ineqc.constraints[1], getfriction(mechanism, ineqc.parentid), mechanism.Δt)
end

function gs(mechanism, ineqc::InequalityConstraint)
    # this is the residual with substracted slacks
    return g(mechanism, ineqc) - ineqc.ssol[2]
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
# function ∂g∂ʳposb(mechanism, ineqc::InequalityConstraint, body::Body)
#     return ∂g∂ʳposb(ineqc.constraints[1], body, nothing)
# end

function ∂g∂ʳvela(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳvela(ineqc.constraints[1], body, nothing, mechanism.Δt)
end
# function ∂g∂ʳvelb(mechanism, ineqc::InequalityConstraint, body::Body)
#     return ∂g∂ʳvelb(ineqc.constraints[1], body, nothing, mechanism.Δt)
# end
