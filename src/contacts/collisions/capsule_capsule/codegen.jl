using Symbolics 
using LinearAlgebra
using SparseArrays
using Scratch

function L_mult(x)
    [x[1] -transpose(x[2:4]); 
     x[2:4] x[1] * I(3) + skew(x[2:4])]
end

# right quaternion multiply as matrix
function R_mult(x)
    [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - skew(x[2:4])]
end

# rotation matrix
function rotation_matrix(q) 
    H = [zeros(1, 3); I(3)]
    transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
end

nz = 4
ny = 6
nw = nz + ny 
nθ = 12

@variables w[1:(nz + ny)] θ[1:nθ] κ[1:1] ϵ[1:1]

function primals(w) 
    return w[1:nz]
end

function duals(w)
    y = w[nz .+ (1:ny)]

    λ = y[1:2] 
    ν = y[2 .+ (1:4)]

    return λ, ν
end

function data(θ)
    pa1 = θ[1:3]
    pa2 = θ[3 .+ (1:3)]
    pb1 = θ[6 .+ (1:3)]
    pb2 = θ[9 .+ (1:3)]

    return pa1, pa2, pb1, pb2
end

function objective(w, θ)
    z = primals(w)
    t = z[1:2] 

    pa1, pa2, pb1, pb2 = data(θ)

    pa = pa1 + t[1] * (pa2 - pa1)
    pb = pb1 + t[2] * (pb2 - pb1) 

    d = pa - pb 

    return dot(d, d)
end

function constraints(w, θ) 
    z = primals(w) 
    t = z[1:2] 
    s = z[2 .+ (1:2)]
    [
        s[1] - (1.0 - t[1]);
        s[2] - (1.0 - t[2]);
    ]
end

function lagrangian(w, θ)
    # initialize
    L = 0.0

    # primals 
    z = primals(w) 

    # duals 
    λ, ν = duals(w)

    # objective 
    J = objective(w, θ) 
    L += J 

    # constraints
    c = constraints(w, θ)
    L += dot(λ, c)

    # inequalities 
    L -= dot(ν, z)

    return L 
end

L = lagrangian(w, θ)
Lz = Symbolics.gradient(L, w[1:nz])

function residual(w, θ, κ)
    # primals 
    z = primals(w) 

    # duals 
    λ, ν = duals(w)

    # Lagrangian 
    lag = lagrangian(w, θ)
    lagz = Symbolics.gradient(lag, z)

    con = constraints(w, θ) 

    comp = ν .* z
           
    res = [
            lagz;
            con;
            comp .- κ;
          ]

    return res 
end

r = residual(w, θ, κ)
rw = Symbolics.jacobian(r, w)
rθ = Symbolics.jacobian(r, θ)

# regularization
for i = 1:nz 
    rw[i, i] += ϵ[1]
end
for j = 1:ny 
    rw[nz + j, nz + j] -= ϵ[1]
end

r_capsule_capsule_func = eval(Symbolics.build_function(r, w, θ, κ)[2])
rw_capsule_capsule_func = eval(Symbolics.build_function(rw, w, θ, ϵ)[2])
rθ_capsule_capsule_func = eval(Symbolics.build_function(rθ, w, θ)[2])

path_collisions = @get_scratch!("collisions")
path_expr = joinpath(path_collisions, "capsule_capsule" * ".jld2")

@save path_expr r_capsule_capsule_func rw_capsule_capsule_func rθ_capsule_capsule_func

# pre-allocate
nθ
r0 = zeros(nw)
rw0 = zeros(nw, nw)
rθ0 = zeros(nw, nθ)

## initialization 
w0 = [0.5; 0.5; 1.0; 1.0; 0.0; 0.0; 1.0; 1.0; 1.0; 1.0]
θ0 = rand(12)
κ0 = [1.0]

r_capsule_capsule_func(r0, w0, θ0, κ0)
rw_capsule_capsule_func(rw0, w0, θ0, [1.0e-5])
rθ_capsule_capsule_func(rθ0, w0, θ0)
cond(rw0)
rank(rw0)
rw0 \ r0
rw0 \ rθ0