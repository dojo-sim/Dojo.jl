"""
$(TYPEDEF)

A `Storage` can be used to store simulation results for a [`Mechanism`](@ref).

Indexing example:
    storage.x[body_id][time_step_number][2] # returns the y-position of the specified body at the specified time step.

# Important attributes
* `x`: Contains the position for each body and each time step.
* `q`: Contains the orientation for each body and each time step.
* `v`: Contains the velocity for each body and each time step.
* `ω`: Contains the angular velocity for each body and each time step.

# Constructors
    Storage{T}(step_range, nbodies)
    Storage(step_range, nbodies)
    Storage{T}(nsteps, nbodies)
    Storage(nsteps, nbodies)
"""
struct Storage{T,N}
    x::Vector{Vector{SVector{3,T}}}
    q::Vector{Vector{UnitQuaternion{T}}}
    v::Vector{Vector{SVector{3,T}}}
    ω::Vector{Vector{SVector{3,T}}}

    function Storage{T}(steps, nbodies) where T
        x = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        q = [[one(UnitQuaternion{T}) for i = steps] for j = 1:nbodies]
        v = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        ω = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        new{T,length(steps)}(x, q, v, ω)
    end

    Storage(steps, nbodies) = Storage{Float64}(steps, nbodies)
    Storage{T}(nsteps::Integer, nbodies) where T = Storage{T}(Base.OneTo(nsteps), nbodies)
    Storage(nsteps::Integer, nbodies) = Storage{Float64}(nsteps, nbodies)

    function Storage(x::Vector{<:Vector{<:AbstractVector{T}}},q::Vector{Vector{UnitQuaternion{T}}}) where T
        steps = Base.OneTo(length(x[1]))
        nbodies = length(x)

        v = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        ω = [[szeros(T, 3) for i = steps] for j = 1:nbodies]

        new{T,length(steps)}(x, q, v, ω)
    end

    Storage{T}() where T = Storage{T}(Base.OneTo(0),0)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, storage::Storage{T,N}) where {T,N}
    summary(io, storage)
end

function downsample(storage::Storage{T,N}, n::Int) where {T,N}
    steps = N ÷ n
    nbodies = length(storage.x)
    s = Storage(steps, nbodies)
    for i = 1:nbodies
        s.x[i] = storage.x[i][1:n:end]
        s.q[i] = storage.q[i][1:n:end]
        s.v[i] = storage.v[i][1:n:end]
        s.ω[i] = storage.ω[i][1:n:end]
    end
    return s
end
