mutable struct Friction{T} <: Component{T}
    id::Int64
    name::String

    Bx::SMatrix{4,3,T,12}
    p::SVector{3,T}

    parentid::Int64 # the body to which the friction constraint is attached
    childids::SVector{2,Int64} # friction and beta bounds

    βsol::Vector{SVector{4,T}}
    γsolref::Vector{SVector{1,T}} #TODO this is a reference to the associated impact γ to avoid allocations

    # B(z) zdot + 1 ψ - η = 0
    d::SVector{4,T} # seems to hold the friction equality

    function Friction(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3), name::String="") where T
        N = 4

        Bx = SA{T}[
            1 0 0
            -1 0 0
            0 1 0
            0 -1 0
        ]

        frictionid = getGlobalID()
        impact = InequalityConstraint(Impact(body, normal; p = p, offset = offset))
        frictionbound = InequalityConstraint(FrictionBound(cf, frictionid, impact.id))
        betabound = InequalityConstraint(BetaBound(frictionid))

        βsol = [szeros(T, N) for i=1:2] # β
        γsolref = impact.γsol # γ

        d = szeros(T, 4)

        new{T}(frictionid, name, Bx, p, body.id, [frictionbound.id; betabound.id], βsol, γsolref, d), [impact; frictionbound; betabound]
    end
end


Base.length(::Friction) = 4

@inline ∂g∂ʳself(mechanism, fric::Friction{T}) where {T} = szeros(T,4,4)

@inline Bq(Bxmat, p, q) = Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p) * 2.0
@inline Bmat(Bxmat, p, q) = [Bxmat Bq(Bxmat, p, q)]


function g(mechanism, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    x, v, q, ω = fullargssol(body.state)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    Bxmat = fric.Bx
    # transforms the velocities of the contact point attached to the link into velocities along all 4 axes of the friction pyramid
    Bqmat = Bq(Bxmat, fric.p, q)

    # this is the B(z) * zdot part of the MDP equation
    # B(z) * zdot + 1 * ψ - η = 0
    fric.d = Bxmat*v + Bqmat*ω

    for childid in fric.childids
        constraintForceMapping!(mechanism, fric, getineqconstraint(mechanism, childid))
    end

    return fric.d
end


# contribution of the inequality constraint (friction) to the dynamics equation d = 0
@inline function constraintForceMapping!(mechanism, body::Body, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    # x1, q1 current state
    x, q = posargsk(body.state)

    # transform the forces βsol expressed along the pyramid axes into the coordinates of the link (x, q)
    body.state.d -= Bmat(fric.Bx, fric.p, q)' * fric.βsol[2]
    return
end

@inline function ∂gab∂ʳba(mechanism, body::Body, fric::Friction)
    x, q = posargsk(body.state)
    B = Bmat(fric.Bx, fric.p, q)
    # jacobian of the friction contribution to the dynamics equation, take wrt β
    # return ∂body.state.d / ∂β; ∂fric.d / ∂zdot
    return -B', B
end


@inline function ∂gab∂ʳba(mechanism, fric::Friction, ineqc::InequalityConstraint)
    G = ∂g∂beta(fric, ineqc.constraints[1])

    return G, -G'
end
