
function divi(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    # adapted from https://www.sciencedirect.com/science/article/pii/S0893965914000469
    # complexity O(n)
    n = length(u)
    @assert length(v) == n
    s = u[2:end] ./ u[1] # n-1
    α = - s' * s # n-1
    t = v - 1 / (1 + α) * (v[1] - s' * v[2:end]) * [0.0; s] # 2n-1
    w = t - [s' * t[2:end]; zeros(n-1)] # n
    w = w ./ u[1] # n
    return w
end

function lift(u::AbstractVector{T}) where {T}
    n = length(u)
    U = Array(Diagonal(u[1] * ones(n)))
    U[2:end,1] = u[2:end]
    U[1,2:end] = u[2:end]
    return U
end

function matinv(u::AbstractVector{T}) where {T}
    n = length(u)
    Ui = zeros(n, n)

    # s = u[2:end] ./ u[1] # n-1
    s = u[2:end] # n-1
    α = - transpose(s) * s # n-1
    S1 = zeros(eltype(u), n, n)
    S2 = zeros(eltype(u), n, n)
    S1[1,2:n] = u[2:end] # s
    S2[2:n,1] = u[2:end] # s
    Ui = (I - S1) * (I - 1/(1+α) * (S2 * (I - S1)))
    return Ui
end

n = 4
u = [1; rand(n-1)]
v = rand(n)
Ui = matinv(u)
norm(Ui * v - lift(u) \ v, Inf)
Ui

@variables u_[1:n], v_[1:n]
u__ = [u_[1], u_[2], u_[3], u_[4]]
Ui_ = matinv(u__)
println(Ui_)




# dims
n = 20
m = 4
p = 5

#data
b = rand(n)
A = [rand(p, p) for i = 1:m]
A = [a'*a for a ∈ A]
B = [rand(p, p) for i = 1:m-1]
C = [rand(p, p) for i = 1:m-1]
M = zeros(n, n)
for i = 1:m
    inds = p*(i-1) .+ (1:p)
    M[inds, inds] = A[i]
    (i < m) && (M[end-p+1:end, inds] = B[i])
    (i < m) && (M[inds, end-p+1:end] = C[i])
end
plot(Gray.(abs.(M)))
x0 = M \ b
x1 = zeros(n)
for i = 1:m-1
    inds = p*(i-1) .+ (1:p)
    b_ = [b[inds]; b[end-p+1:end] ./ 3]
    M_ = [A[i] C[i]; B[i] A[m] ./ 3]
    x_ = M_ \ b_
    x1[inds] = x_[1:p]
    x1[end-p+1:end] += x_[p+1:2p]
end
norm(x0 - x1)
x0
x1
M * x0 - b
M * x1 - b







function res(x3, q3)
    Δt = 0.01
    ainv3 = [0 0 1.]
    p = [0.25, 0.25, -0.25]
    offset = [0, 0, 0.]
    return ainv3 * (x3 + vrotate(p,q3) - offset)
end

function dresdx3(x3, q3)
    Δt = 0.01
    ainv3 = [0 0 1.]
    p = [0.25, 0.25, -0.25]
    offset = [0, 0, 0.]
    return ainv3
end

function dresdv2(x2, v2, q2, ω2)
    Δt = 0.01
    ainv3 = [0 0 1.]
    p = [0.25, 0.25, -0.25]
    offset = [0, 0, 0.]

    x3 = x2 + v2 * Δt
    q3 = q2 * ωbar(ω2, Δt) * Δt / 2
    return dresdx3(x3, q3) * Δt
end

function dresdq3(x3, q3)
    Δt = 0.01
    ainv3 = [0 0 1.]
    p = [0.25, 0.25, -0.25]
    offset = [0, 0, 0.]
    # return (ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(p)))) #* LVᵀmat(q3)
    return ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(p)))
end

function dresdω2(x2, v2, q2, ω2)
    Δt = 0.01
    ainv3 = [0 0 1.]
    p = [0.25, 0.25, -0.25]
    offset = [0, 0, 0.]

    x3 = x2 + v2 * Δt
    q3 = q2 * ωbar(ω2, Δt) * Δt / 2
    return dresdq3(x3, q3) * Lmat(q2) * derivωbar(SVector{3}(ω2[1], ω2[2], ω2[3]), Δt) * Δt / 2
end

function res(x2, v2, q2, ω2)
    x3 = x2 + v2 * Δt
    q3 = q2 * ωbar(ω2, Δt) * Δt / 2
    return res(x3, q3)
end

Random.seed!(10)
x3 = rand(3)
q3 = UnitQuaternion(rand(3))
q3v = [q3.w, q3.x, q3.y, q3.z]

r30 = res(x3, q3)
drdx30 = dresdx3(x3, q3)
drdq30 = dresdq3(x3, q3)

drdx31 = fdjac(x -> res(x, q3), x3)
drdq31 = fdjac(q -> res(x3, UnitQuaternion(q, false)), q3v)
norm(drdx30 - drdx31)
# norm(drdq30 - drdq31 * LVᵀmat(q3))
norm(drdq30 - drdq31)

Random.seed!(10)
x2 = rand(3)
v2 = rand(3)
q2 = UnitQuaternion(rand(3))
q2v = [q2.w, q2.x, q2.y, q2.z]
ω2 = rand(3)

r20 = res(x2, v2, q2, ω2)
drdv20 = dresdv2(x2, v2, q2, ω2)
drdω20 = dresdω2(x2, v2, q2, ω2)

drdv21 = fdjac(v -> res(x2, v, q2, ω2), v2)
drdω21 = fdjac(ω -> res(x2, v2, q2, ω), ω2)

norm(drdv20 - drdv21)
norm(drdω20 - drdω21)


function fdjac(f, x; δ = 1e-5)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end
















ineqc = mech.ineqconstraints[3]
cont = ineqc.constraints[1]
p = cont.p
Bx = cont.Bx

function Bv(q2, ω2)
    Bx * VRᵀmat(q2) * LVᵀmat(q2) * 2 * skew(-p) * ω2
end

Random.seed!(1)
q2 = UnitQuaternion(rand(4)...)
ω2 = rand(3)

# q2 = storage.q[1][end]
# ω2 = storage.ω[1][end]
Bv(q2, ω2)

function fdjac(f, x; δ = 1e-5)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end

q2v = [q2.w, q2.x, q2.y, q2.z]
fdjac(q -> Bv(UnitQuaternion(q, false), ω2), q2v) * G(q2v)
dBqω(q2v, ω2, cont.p)





function link_sizes(vertices::Vector{Int}, links::Vector{V}, vsizes::Vector{Int}) where {V}
    lsizes = Vector{Vector{Int}}()
    for l in links
        s1 = vsizes[vertices[findfirst(x -> x == l[1], vertices)]]
        s2 = vsizes[vertices[findfirst(x -> x == l[2], vertices)]]
        push!(lsizes, [s1, s2])
    end
    return lsizes
end

function complete_matrix(vertices, links, vsizes, lsizes, Mvertices, Mlinks)
    nv = length(vertices)
    nl = length(links)
    n = sum(vsizes)
    M = zeros(n, n)

    off = 0
    for i = 1:nv
        ni = vsizes[i]
        M[off .+ (1:ni), off .+ (1:ni)] = Mvertices[i]
        off += ni
    end
    off = 0
    for i = 1:nl
        ni = vsizes[i]
        n1, n2 = lsizes[i]
        @show n1, n2
        M[off + ni .+ (1:n2), off .+ (1:n1)] = Mlinks[i][1]
        M[off .+ (1:n1), off + ni .+ (1:n2)] = Mlinks[i][2]
        off += ni
    end
    return M
end



vertices = [1, 2, 3, 4]
links = [[1,2], [2,3], [3,4]]
vsizes = [3, 5, 6, 2]
lsizes = link_sizes(vertices, links, vsizes)
Mvertices = [0.0*rand(s,s) for s in vsizes]
Mlinks = [[rand(s[2], s[1]), rand(s[1], s[2])] for s in lsizes]
M = complete_matrix(vertices, links, vsizes, lsizes, Mvertices, Mlinks)
plot(Gray.(1e10 .* abs.(M)))

Mlinks
