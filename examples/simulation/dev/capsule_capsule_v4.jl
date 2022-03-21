using Symbolics 
using RoboDojo
using LinearAlgebra
using SparseArrays

using MeshCat
using GeometryBasics
using Colors
using CoordinateTransformations
using Rotations

# visualize
function set_background!(vis::Visualizer; top_color=RGBA(1,1,1.0), bottom_color=RGBA(1,1,1.0))
    setprop!(vis["/Background"], "top_color", top_color)
    setprop!(vis["/Background"], "bottom_color", bottom_color)
end
vis = Visualizer()
open(vis)
set_background!(vis)

function hat(x)
    return [0 -x[3] x[2];
            x[3] 0 -x[1];
           -x[2] x[1] 0]
end

function L_mult(x)
    [x[1] -transpose(x[2:4]); 
     x[2:4] x[1] * I(3) + hat(x[2:4])]
end

# right quaternion multiply as matrix
function R_mult(x)
    [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - hat(x[2:4])]
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
@variables z[1:nz] y[1:ny] θ[1:nθ] κ[1:1]
w = [z; y]

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

    return dot(d, d) + 1.0e-5 * (t[1] - 0.5)^2 + 1.0e-5 * (t[2] - 0.5)^2
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
Lz = Symbolics.gradient(L, z)

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

r_func = eval(Symbolics.build_function(r, w, θ, κ)[2])
rw_func = eval(Symbolics.build_function(rw, w, θ)[2])
rθ_func = eval(Symbolics.build_function(rθ, w, θ)[2])

# pre-allocate
r0 = zeros(nw)
rw0 = zeros(nw, nw)
rθ0 = zeros(nw, nθ)

function rw_func_reg(rw, w, θ) 
    rw_func(rw, w, θ)
    rw .+= Diagonal([1.0e-5 * ones(nz); -1.0e-5 * ones(ny)])
    return 
end


## setup 
xa = [-1.0; 0.0; 0.0]
qa = Quaternion(RotY(0.0 * π) * RotX(0.0 * π))
qa = [qa.w; qa.x; qa.y; qa.z]
ra = 0.1
ha = 0.2

pa1 = xa + rotation_matrix(qa) * [0.0; 0.0;  0.5 * ha]
pa2 = xa + rotation_matrix(qa) * [0.0; 0.0; -0.5 * ha]

xb = [1.0; 0.0; 0.0]
qb = Quaternion(RotZ(0.0 * π) * RotY(0.0 * π) * RotX(0.0))
qb = [qb.w; qb.x; qb.y; qb.z]
rb = 0.1 
hb = 0.2
pb1 = xb + rotation_matrix(qb) * [0.0; 0.0;  0.5 * hb]
pb2 = xb + rotation_matrix(qb) * [0.0; 0.0; -0.5 * hb]

## initialization 
w0 = [0.5; 0.5; 1.0; 1.0; 0.0; 0.0; 1.0; 1.0; 1.0; 1.0]
θ0 = [pa1; pa2; pb1; pb2]
κ0 = [1.0]

r_func(r0, w0, θ0, κ0)
rw_func_reg(rw0, w0, θ0)
rθ_func(rθ0, w0, θ0)
cond(rw0)
rank(rw0)
rw0 \ r0
rw0 \ rθ0

## solver 
idx = RoboDojo.IndicesOptimization(
    nw, nw, 
    [collect(1:4), collect(7:10)], [collect(1:4), collect(7:10)],
    Vector{Vector{Vector{Int}}}(), Vector{Vector{Vector{Int}}}(),
    collect(1:6),
    collect(7:10),
    Vector{Int}(),
    Vector{Vector{Int}}(),
    collect(7:10),
)

ip = RoboDojo.interior_point(w0, θ0;
    s = RoboDojo.Euclidean(length(w0)),
    idx = idx,
    r! = r_func, 
    rz! = rw_func_reg, 
    rθ! = rθ_func,
    r  = zeros(idx.nΔ),
    rz = zeros(idx.nΔ, idx.nΔ),
    rθ = zeros(idx.nΔ, length(θ0)),
    opts = RoboDojo.InteriorPointOptions(
            undercut=5.0,
            γ_reg=0.1,
            r_tol=1e-5,
            κ_tol=1e-5,  
            max_ls=50,
            ϵ_min=0.05,
            diff_sol=true,
            verbose=true))

RoboDojo.interior_point_solve!(ip)

ta = ip.z[1]
tb = ip.z[2]
pa = pa1 + ta * (pa2 - pa1)
pb = pb1 + tb * (pb2 - pb1)


# visualize 
vis = Visualizer() 
open(vis) 

# capsule a
cyl1 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * ha), Point(0.0, 0.0, 0.5 * ha), ra)
sph1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), ra)
color1 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap1][:cyl], cyl1, MeshPhongMaterial(color=color1))
setobject!(vis[:cap1][:base], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:base], Translation(0.0, 0.0, -0.5 * ha))
setobject!(vis[:cap1][:tip], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:tip], Translation(0.0, 0.0, 0.5 * ha))

# capsule b
cyl2 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * hb), Point(0.0, 0.0, 0.5 * hb), rb)
sph2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), rb)
color2 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap2][:cyl], cyl2, MeshPhongMaterial(color=color2))
setobject!(vis[:cap2][:base], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:base], Translation(0.0, 0.0, -0.5 * hb))
setobject!(vis[:cap2][:tip], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:tip], Translation(0.0, 0.0, 0.5 * hb))

# set configuration
settransform!(vis[:cap1], compose(Translation(xa), LinearMap(Quaternion(qa..., false))))
settransform!(vis[:cap2], compose(Translation(xb), LinearMap(Quaternion(qb..., false))))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(pa - dir * 0.0 * ra))
settransform!(vis[:cs2], Translation(pb + dir * 0.0 * rb))