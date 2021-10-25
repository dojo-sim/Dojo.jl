"""
$(TYPEDEF)

A `Mechanism` contains the [`Origin`](@ref), [`Body`](@ref)s, and [`EqualityConstraint`](@ref)s of a system and can be used for simulation.
# Important attributes
* `origin`:        The origin of a mechanism.
* `bodies`:        The bodies of a mechanism (Dict).
* `eqconstraints`: The equality constraints (joints) of a mechanism (Dict).
* `Δt`:            The time step of the mechanism.
* `g`:             The gravitational acceleration in z-direction.

# Constuctors
    Mechanism(origin, bodies; Δt, g)
    Mechanism(origin, bodies, eqcs; Δt, g)
    Mechanism(urdf_filename; floating, Δt, g)
"""
mutable struct Mechanism{T,Nn,Ne,Nb,Nf,Ni} <: AbstractMechanism{T,Nn,Ne,Nb,Nf,Ni}
    origin::Origin{T}
    eqconstraints::UnitDict{Base.OneTo{Int64},<:EqualityConstraint{T}}
    bodies::UnitDict{UnitRange{Int64},Body{T}}
    frictions::UnitDict{UnitRange{Int64},<:Friction{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    system::System{Nn}
    residual_entries::Vector{Entry}
    matrix_entries::SparseMatrixCSC{Entry, Int64}
    diagonal_inverses::Vector{Entry}

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T
    rvio::T
    bvio::T
    ν::T
    νaff::T

    Δt::T
    g::T

    α::T
    μ::T


    function Mechanism(origin::Origin{T},bodies::Vector{<:Body{T}},
            eqcs::Vector{<:EqualityConstraint{T}}, ineqcs::Vector{<:InequalityConstraint{T}},
            frics::Vector{<:Friction{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        resetGlobalID()
        order = getGlobalOrder()

        for body in bodies
            initknotpoints!(body.state, order)
            if norm(body.m) == 0 || norm(body.J) == 0
                @info "Potentially bad inertial properties detected"
            end
        end

        Ne = length(eqcs)
        Nb = length(bodies)
        Nf = length(frics)
        Ni = length(ineqcs)
        Nn = Ne + Nb + Nf + Ni

        nodes = [eqcs;bodies;frics;ineqcs]
        oldnewid = Dict([node.id=>i for (i,node) in enumerate(nodes)]...)

        for node in nodes
            node.id = oldnewid[node.id]
            # @show node.id
            # @show typeof(node).name
            # error("gere")
            if typeof(node) <: Union{AbstractConstraint, Friction}
                # @show typeof(node).name
                # @show node.parentid
                # @show node.childids
                node.parentid = get(oldnewid, node.parentid, nothing)
                node.childids = [get(oldnewid, childid, nothing) for childid in node.childids]
                # @show node.parentid
                # @show node.childids
            end
        end

        system = create_system(origin, eqcs, bodies, frics, ineqcs)
        residual_entries = deepcopy(system.vector_entries)
        matrix_entries = deepcopy(system.matrix_entries)
        diagonal_inverses = deepcopy(system.diagonal_inverses)

        normf = 0
        normΔs = 0
        rvio = 0
        bvio = 0
        ν = 0
        νaff = 0

        eqcs = UnitDict(eqcs)
        bodies = UnitDict((bodies[1].id):(bodies[Nb].id), bodies)
        Nf > 0 ? (frics = UnitDict((frics[1].id):(frics[Nf].id), frics)) : (frics = UnitDict(0:0, frics))
        Ni > 0 ? (ineqcs = UnitDict((ineqcs[1].id):(ineqcs[Ni].id), ineqcs)) : (ineqcs = UnitDict(0:0, ineqcs))


        α = 1
        μ = 1

        new{T,Nn,Ne,Nb,Nf,Ni}(origin, eqcs, bodies, frics, ineqcs, system, residual_entries, matrix_entries, diagonal_inverses,
            normf, normΔs, rvio, bvio, ν, νaff, Δt, g, α, μ)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{<:Body{T}},eqcs::Vector{<:EqualityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        ineqcs = InequalityConstraint{T}[]
        frics = Friction{T}[]
        return Mechanism(origin, bodies, eqcs, ineqcs, frics; Δt = Δt, g = g)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{<:Body{T}},eqcs::Vector{<:EqualityConstraint{T}},ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        frics = Friction{T}[]
        return Mechanism(origin, bodies, eqcs, ineqcs, frics; Δt = Δt, g = g)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{<:Body{T}},ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, ineqcs; Δt = Δt, g = g)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{<:Body{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, Δt = Δt, g = g)
    end

    function Mechanism(filename::AbstractString; floating::Bool=false, type::Type{T} = Float64, Δt::Real = .01, g::Real = -9.81) where T
        origin, links, joints, loopjoints = parse_urdf(filename, floating, T)
        mechanism = Mechanism(origin, links, [joints;loopjoints], Δt = Δt, g = g)
        set_parsed_values!(mechanism, loopjoints)

        return mechanism
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, mechanism::AbstractMechanism{T,Nn,Ne,Nb,Nf,0}) where {T,Nn,Ne,Nb,Nf}
    summary(io, mechanism)
    println(io, " with ", Nb, " bodies and ", Ne, " constraints")
    println(io, " Δt: "*string(mechanism.Δt))
    println(io, " g:  "*string(mechanism.g))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, mechanism::AbstractMechanism{T,Nn,Ne,Nb,Nf,Ni}) where {T,Nn,Ne,Nb,Nf,Ni}
    summary(io, mechanism)
    println(io, " with ", Nb, " bodies, ", Ne, " equality constraints, and ", Ni, " inequality constraints")
    println(io, " Δt: "*string(mechanism.Δt))
    println(io, " g:  "*string(mechanism.g))
end
