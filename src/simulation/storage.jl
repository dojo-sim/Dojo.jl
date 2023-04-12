"""
    Storage{T,N}

    contains maximal-representation trajectories

    x: position 
    q: orientation (Quaternion)
    v: linear velocity (midpoint) 
    ω: angular velocity (midpoint)
    px: linear momentum
    pq: angular momentum 
    vl: linear velocity
    ωl: angular velocity
"""
struct Storage{T,N}
    x::Vector{Vector{SVector{3,T}}}
    q::Vector{Vector{Quaternion{T}}}
    v::Vector{Vector{SVector{3,T}}}
    ω::Vector{Vector{SVector{3,T}}}
    px::Vector{Vector{SVector{3,T}}}
    pq::Vector{Vector{SVector{3,T}}}
    vl::Vector{Vector{SVector{3,T}}}
    ωl::Vector{Vector{SVector{3,T}}}

    function Storage{T}(steps, nbodies) where T
        x = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        q = [[one(Quaternion{T}) for i = steps] for j = 1:nbodies]
        v = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        ω = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        px = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        pq = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        vl = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        ωl = [[szeros(T, 3) for i = steps] for j = 1:nbodies]
        new{T,length(steps)}(x, q, v, ω, px, pq, vl, ωl)
    end

    Storage(steps, nbodies) = Storage{Float64}(steps, nbodies)
    Storage{T}(nsteps::Integer, nbodies) where T = Storage{T}(Base.OneTo(nsteps), nbodies)
    Storage(nsteps::Integer, nbodies) = Storage{Float64}(nsteps, nbodies)

    Storage{T}() where T = Storage{T}(Base.OneTo(0),0)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, storage::Storage{T,N}) where {T,N}
    println(io, "Storage for "*string(N)*" steps of "*string(length(storage.x))*" bodies.")
end

Base.length(storage::Storage{T,N}) where {T,N} = N

function save_to_storage!(mechanism::Mechanism, storage::Storage, i::Int)
    for (ind, body) in enumerate(mechanism.bodies)
        state = body.state
        storage.x[ind][i] = state.x2 # x2
        storage.q[ind][i] = state.q2 # q2
        storage.v[ind][i] = state.v15 # v1.5
        storage.ω[ind][i] = state.ω15 # ω1.5
        q2 = state.q2
        px2, pq2 = momentum(mechanism, body) # p in world frame
        v2 = px2 / body.mass # in world frame
        ω2 = body.inertia \ (vector_rotate(pq2, inv(q2))) # in body frame, we rotate using the current quaternion q2 = state.q2
        storage.px[ind][i] = px2 # px2
        storage.pq[ind][i] = pq2 # pq2
        storage.vl[ind][i] = v2 # v2
        storage.ωl[ind][i] = ω2 # ω2
    end
    return
end

function generate_storage(mechanism::Mechanism, z)
    N = length(z)
    M = length(mechanism.bodies)
	storage = Storage{Float64}(N, M)

    for t = 1:N
        off = 0
        for i = 1:M 
            storage.x[i][t] = z[t][off .+ (1:3)]
            storage.v[i][t] = z[t][off .+ (4:6)]
            storage.q[i][t] = Quaternion(z[t][off .+ (7:10)]...)
            storage.ω[i][t] = z[t][off .+ (11:13)]
            off += 13
        end
    end

    return storage
end

function get_maximal_state(storage::Storage{T,N}) where {T,N}
	z = [get_maximal_state(storage, i) for i=1:N]
	return z
end

function get_maximal_state(storage::Storage{T,N}, i::Int) where {T,N}
	Nb = length(storage.x)
	z = zeros(13 * Nb)
	for j = 1:Nb
		x2 = storage.x[j][i]
		q2 = storage.q[j][i]
		v15 = storage.v[j][i]
		ω15 = storage.ω[j][i]
		z[13 * (j-1) .+ (1:13)] = [x2; v15; vector(q2); ω15]
	end
	return z
end
