mutable struct EqualityConstraint{T,N,Nc,Cs} <: Constraint{T,N}
    id::Int64
    name::String
    isspring::Bool
    isdamper::Bool

    constraints::Cs
    parentid::Union{Int64,Nothing}
    childids::SVector{Nc,Int64}
    inds::SVector{Nc,SVector{2,Int64}} # indices for minimal coordinates, assumes joints # Nc = 2 THIS IS SPECIAL CASED

    λsol::Vector{SVector{N,T}}

    function EqualityConstraint(data; name::String="")
        jointdata = Tuple{Joint,Int64,Int64}[]
        for info in data
            push!(jointdata, info)
        end

        T = getT(jointdata[1][1])

        isspring = false
        isdamper = false
        parentid = jointdata[1][2]
        childids = Int64[]
        constraints = Joint{T}[]
        inds = Vector{Int64}[]
        N = 0
        for set in jointdata
            set[1].spring != 0 && (isspring = true)
            set[1].damper != 0 && (isdamper = true)

            push!(constraints, set[1])
            @assert set[2] == parentid
            push!(childids, set[3])

            Nset = ηlength(set[1])
            if isempty(inds)
                push!(inds, [1;3-Nset])
            else
                push!(inds, [last(inds)[2]+1; last(inds)[2]+3-Nset])
            end
            N += Nset
        end
        constraints = Tuple(constraints)
        Nc = length(constraints)
        λsol = [zeros(T, N) for i=1:2]
        return new{T,N,Nc,typeof(constraints)}(getGlobalID(), name, isspring, isdamper, constraints, parentid, childids, inds, λsol)
    end
end

# joint_limits_length(eqc::EqualityConstraint) = sum(joint_limits_length.(eqc.constraints))

function setPosition!(mechanism, eqc::EqualityConstraint, xθ; iter::Bool=true)
    if !iter
        _setPosition!(mechanism, eqc, xθ)
    else
        currentvals = minimalCoordinates(mechanism)
        _setPosition!(mechanism, eqc, xθ)
        for id in recursivedirectchildren!(mechanism.system, eqc.id)
            component = getcomponent(mechanism, id)
            if component isa EqualityConstraint
                _setPosition!(mechanism, component, currentvals[id])
            end
        end
    end

    return
end

# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
function _setPosition!(mechanism, eqc::EqualityConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    @assert length(xθ)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.parentid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δx = getPositionDelta(eqc.constraints[i], body1, body2, xθ[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])]) # in body1's frame
        Δq = getPositionDelta(eqc.constraints[i+1], body1, body2, xθ[SUnitRange(eqc.inds[i+1][1], eqc.inds[i+1][2])]) # in body1's frame

        p1, p2 = eqc.constraints[i].vertices
        setPosition!(body1, body2; p1=p1, p2=p2, Δx=Δx, Δq=Δq)
    end
    return
end

function setVelocity!(mechanism, eqc::EqualityConstraint{T,N,Nc}, vω) where {T,N,Nc}
    # vω is already in body1 frame
    @assert length(vω)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.parentid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δv = getVelocityDelta(eqc.constraints[i], body1, body2, vω[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])]) # projection in body1 frame
        Δω = getVelocityDelta(eqc.constraints[i+1], body1, body2, vω[SUnitRange(eqc.inds[i+1][1], eqc.inds[i+1][2])]) # projection in body1 frame
        p1, p2 = eqc.constraints[i].vertices
        setVelocity!(body1, body2; p1=p1, p2=p2, Δv=Δv, Δω=Δω)
    end
    return
end

function setForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ::AbstractVector) where {T,N,Nc}
    @assert length(Fτ)==controldim(eqc)
    for i = 1:Nc
        r_idx = SUnitRange(eqc.inds[i][1], eqc.inds[i][2])
        length(r_idx) == 0 && continue
        setForce!(eqc.constraints[i], Fτ[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])])
    end
    return
end

function addForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ::AbstractVector) where {T,N,Nc}
    @assert length(Fτ)==controldim(eqc)
    for i = 1:Nc
        addForce!(eqc.constraints[i], Fτ[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])])
    end
    return
end

"""
    minimalCoordinates(mechanism, eqconstraint)

Gets the minimal coordinates of joint `eqconstraint`.
"""
@generated function minimalCoordinates(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimalCoordinates(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@generated function minimalVelocities(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimalVelocities(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function constraintForceMapping!(mechanism, body::Body, eqc::EqualityConstraint)
    # @show size(∂g∂ʳpos(mechanism, eqc, body))
    # @show size(eqc.λsol[2])
    # @warn "we have an issue here we need to filter out the si's that are inside λsol because they do not factor in the dynamics"
    # λ = []
    # for (i, joint) in enumerate(eqc.constraints)
    #     # Nli = joint_limits_length(joint)
    #     λi = eqc.λsol[2][λindex(eqc, i)]
    #     si, γi = get_sγ(joint, λi)
    #     push!(λ, λi[1:length(joint)])
    #     push!(λ, γi)
    #     # push!(comp_jac, [zeros(2Nli, length(joint)) Diagonal(γi) Diagonal(si)])
    # end

    body.state.d -= zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2]

    # if length(eqc.λsol[2]) == 3 && body.id ∈ eqc.childids
    #     Fτ = zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body)) * eqc.λsol[2]
    # end
    eqc.isspring && (body.state.d -= springforce(mechanism, eqc, body))
    eqc.isdamper && (body.state.d -= damperforce(mechanism, eqc, body))
    return
end

@inline function ∂constraintForceMapping!(mechanism, body::Body, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    if body.id == eqc.parentid
        _dGa!(mechanism, body, eqc)
    elseif body.id ∈ eqc.childids
        _dGb!(mechanism, body, eqc)
    else
        error()
    end
    return
end

function _dGa!(mechanism, pbody::Body, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    Δt = mechanism.Δt
    _, _, q2, ω2 = fullargssol(pbody.state)
    M = ∂integration(q2, ω2, Δt)

    off = 0
    for i in 1:Nc
        joint = eqc.constraints[i]
        Aᵀ = zerodimstaticadjoint(constraintmat(joint))
        Nj = length(joint)
        cbody = getbody(mechanism, eqc.childids[i])
        eqc.isspring && (pbody.state.D -= ∂springforcea∂vela(joint, pbody, cbody, Δt))
        eqc.isdamper && (pbody.state.D -= ∂damperforcea∂vela(joint, pbody, cbody, Δt))
        off += Nj
    end
    return nothing
end

function _dGb!(mechanism, cbody::Body, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    Δt = mechanism.Δt
    x2, v2, q2, ω2 = fullargssol(cbody.state)
    M = ∂integration(q2, ω2, Δt)

    off = 0
    for i in 1:Nc
        if eqc.childids[i] == cbody.id
            joint = eqc.constraints[i]
            Aᵀ = zerodimstaticadjoint(constraintmat(joint))
            Nj = length(joint)
            pbody = getbody(mechanism, eqc.parentid)
            eqc.isspring && (cbody.state.D -= ∂springforceb∂velb(joint, pbody, cbody, Δt))
            eqc.isdamper && (cbody.state.D -= ∂damperforceb∂velb(joint, pbody, cbody, Δt))
        end
    end
    return nothing
end

@generated function g(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]), eqc.λsol[2][λindex(eqc,$i)], mechanism.Δt)) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function springforce(mechanism, eqc::EqualityConstraint, body::Body)
    body.id == eqc.parentid ? (return springforcea(mechanism, eqc, body)) : (return springforceb(mechanism, eqc, body))
end

@inline function damperforce(mechanism, eqc::EqualityConstraint, body::Body)
    body.id == eqc.parentid ? (return damperforcea(mechanism, eqc, body)) : (return damperforceb(mechanism, eqc, body))
end

@inline function ∂gab∂ʳba(mechanism, body::Body{T}, eqc::EqualityConstraint{T,N}) where {T,N}
    # Nl = joint_limits_length(eqc)
    # Z = szeros(T,2Nl,6)
    # @warn "this is wrong we need to put zeros for the columns correcponding to the si's"
    # ∇pos = zeros(T,N,6)
    # ∇vel = zeros(T,N,6)
    # @show "ffff"
    return -∂g∂ʳpos(mechanism, eqc, body)', ∂g∂ʳvel(mechanism, eqc, body)
    # return ∇pos', ∇vel
end

@inline function ∂gab∂ʳba(mechanism, eqc::EqualityConstraint{T,N}, body::Body{T}) where {T,N}
    # Nl = joint_limits_length(eqc)
    # Z = szeros(T,2Nl,6)
    # ∇pos = zeros(T,N,6)
    # ∇vel = zeros(T,N,6)
    # @show "gggg"
    # @warn "this is wrong we need to put zeros for the columsn correcponding to the si's"
    return ∂g∂ʳvel(mechanism, eqc, body), -∂g∂ʳpos(mechanism, eqc, body)'
    # return ∇vel, ∇pos'
end
# @inline ∂gab∂ʳba(mechanism, body::Body, eqc::EqualityConstraint) = -∂g∂ʳpos(mechanism, eqc, body)', ∂g∂ʳvel(mechanism, eqc, body)
# @inline ∂gab∂ʳba(mechanism, eqc::EqualityConstraint, body::Body) = ∂g∂ʳvel(mechanism, eqc, body), -∂g∂ʳpos(mechanism, eqc, body)'


@generated function ∂g∂ʳposa(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂g∂ʳposa(eqc.constraints[$i], body, getbody(mechanism, eqc.childids[$i]), eqc.childids[$i], eqc.λsol[2][λindex(eqc,$i)], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end
@generated function ∂g∂ʳposb(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂g∂ʳposb(eqc.constraints[$i], getbody(mechanism, eqc.parentid), body, eqc.childids[$i], eqc.λsol[2][λindex(eqc,$i)], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂ʳself(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(∂g∂ʳself(eqc.constraints[$i], eqc.λsol[2])) for i = 1:Nc]
    return :(cat($(vec...), dims=(1,2)))
end

@generated function ∂g∂ʳvela(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂g∂ʳvela(eqc.constraints[$i], body, getbody(mechanism, eqc.childids[$i]), eqc.childids[$i], eqc.λsol[2][λindex(eqc,$i)], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂ʳvelb(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂g∂ʳvelb(eqc.constraints[$i], getbody(mechanism, eqc.parentid), body, eqc.childids[$i], eqc.λsol[2][λindex(eqc,$i)], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@inline function springforcea(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += springforcea(eqc.constraints[i], body, getbody(mechanism, eqc.childids[i]), mechanism.Δt, eqc.childids[i])
    end
    return vec
end

@inline function springforceb(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += springforceb(eqc.constraints[i], getbody(mechanism, eqc.parentid), body, mechanism.Δt, eqc.childids[i])
    end
    return vec
end

@inline function damperforcea(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damperforcea(eqc.constraints[i], body, getbody(mechanism, eqc.childids[i]), mechanism.Δt, eqc.childids[i])
    end
    return vec
end
@inline function damperforceb(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damperforceb(eqc.constraints[i], getbody(mechanism, eqc.parentid), body, mechanism.Δt, eqc.childids[i])
    end
    return vec
end

@generated function ∂Fτ∂ua(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂Fτ∂ua(eqc.constraints[$i], body, getbody(mechanism, eqc.childids[$i]), mechanism.Δt, eqc.childids[$i])) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@generated function ∂Fτ∂ub(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(∂Fτ∂ub(eqc.constraints[$i], getbody(mechanism, eqc.parentid), body, mechanism.Δt, eqc.childids[$i])) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@inline function applyFτ!(eqc::EqualityConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    for i=1:Nc
        applyFτ!(eqc.constraints[i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[i]), mechanism.Δt, clear)
    end
    return
end

function Base.cat(eqc1::EqualityConstraint{T,N1,Nc1}, eqc2::EqualityConstraint{T,N2,Nc2}) where {T,N1,N2,Nc1,Nc2}
    @assert eqc1.parentid == eqc2.parentid "Can only concatenate constraints with the same parentid"
    parentid = eqc1.parentid
    if parentid === nothing
        parentid = -1
        nothingflag = true
    else
        nothingflag = false
    end

    constraints = [[eqc1.constraints[i] for i=1:Nc1]; [eqc2.constraints[i] for i=1:Nc2]]
    childids = [[eqc1.childids[i] for i=1:Nc1]; [eqc2.childids[i] for i=1:Nc2]]

    eqc = EqualityConstraint([(constraints[i],parentid,childids[i]) for i=1:Nc1+Nc2]..., name="combined_"*eqc1.name*"_and_"*eqc2.name)
    nothingflag && (eqc.parentid = nothing)

    return eqc
end

function set_spring_damper!(eqcs, spring, damper)
    i = 1
    for eqc in eqcs
        eqc.parentid === nothing && continue
        k = (length(spring) > 1) ? spring[i] : spring
        b = (length(damper) > 1) ? damper[i] : damper
        eqc.isspring = k > 0.0
        eqc.isdamper = b > 0.0
        for joint in eqc.constraints
            joint.spring = max(0.0, k)
            joint.damper = max(0.0, b)
        end
        i += 1
    end
    return eqcs
end

@inline function ∂g∂ʳpos(mechanism, constraint::Constraint, body::Body)
    if body.id == constraint.parentid
        return ∂g∂ʳposa(mechanism, constraint, body)
    else
        return ∂g∂ʳposb(mechanism, constraint, body)
    end
end

@inline function ∂g∂ʳvel(mechanism, constraint::Constraint, body::Body)
    body.id == constraint.parentid ? (return ∂g∂ʳvela(mechanism, constraint, body)) : (return ∂g∂ʳvelb(mechanism, constraint, body))
end

function λindex(eqc::EqualityConstraint{T,N,Nc,Cs}, i::Int) where {T,N,Nc,Cs}
    i0 = 1
    i1 = 0
    for j = 1:i
        i0 = i1 + 1
        joint = eqc.constraints[j]
        i1 += λlength(joint)
    end

    ind = SVector{i1-i0+1,Int}(i0:i1...)
    return ind
end

function resetVars!(eqc::EqualityConstraint{T,N,Nc,Cs}; scale::T=1.0) where {T,N,Nc,Cs}
    λ = []
    for (i, joint) in enumerate(eqc.constraints)
        Nλ = length(joint)
        Nb = blength(joint)
        push!(λ, [szeros(Nλ); scale * sones(2Nb)])
    end
    @show size(vcat(λ...))
    @show size(eqc.λsol[1])
    eqc.λsol[1] = vcat(λ...)
    eqc.λsol[2] = vcat(λ...)
    return
end
