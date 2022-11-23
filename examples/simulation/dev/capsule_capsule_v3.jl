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

nz = 10
ny = 12
nw = nz + ny 
nθ = 12
@variables z[1:nz] y[1:ny] θ[1:nθ] κ[1:1]

function unpack_primals(z) 
    pa = z[0 .+ (1:3)]
    pb = z[3 .+ (1:3)]
    ta = z[6 .+ (1:1)]
    tb = z[7 .+ (1:1)]
    sa = z[8 .+ (1:1)]
    sb = z[9 .+ (1:1)]
    return pa, pb, ta, tb, sa, sb
end

function unpack_duals(y)
   λa = y[0 .+ (1:4)]
   λb = y[4 .+ (1:4)]

   νa = y[8 .+ (1:2)]
   νb = y[10 .+ (1:2)]

   return λa, λb, νa, νb
end

function unpack_data(θ)
    pa1 = θ[0 .+ (1:3)]
    pa2 = θ[3 .+ (1:3)]

    pb1 = θ[6 .+ (1:3)]
    pb2 = θ[9 .+ (1:3)]
    
    return pa1, pa2, pb1, pb2
end

function objective(z, θ)
    pa, pb, ta, tb, sa, sb = unpack_primals(z)
    d = pa - pb
    return dot(d, d) + 1.0e-3 * (ta[1] - 0.5)^2 + 1.0e-3 * (tb[1] - 0.5)^2
end

function constraints(z, θ) 
    pa, pb, ta, tb, sa, sb = unpack_primals(z)
    pa1, pa2, pb1, pb2 = unpack_data(θ)

    [
        pa - (pa2 + ta[1] * (pa1 - pa2));
        sa[1] - (1.0 - ta[1]);
        pb - (pb2 + tb[1] * (pb1 - pb2));
        sb[1] - (1.0 - tb[1]);
    ]
end

function lagrangian(z, y, θ)
    # initialize
    L = 0.0

    # primals 
    pa, pb, ta, tb, sa, sb = unpack_primals(z)

    # duals 
    λa, λb, νa, νb = unpack_duals(y)

    # objective 
    J = objective(z, θ) 
    L += J 

    # constraints
    c = constraints(z, θ)
    L += sum([λa[1] λa[2] λa[3] λa[4] λb[1] λb[2] λb[3] λb[4]] * c)

    # inequalities 
    L -= ta[1] * νa[1]
    L -= sa[1] * νa[2]

    L -= tb[1] * νb[1] 
    L -= sb[1] * νb[2] 

    return L 
end

L = lagrangian(z, y, θ)
Lz = Symbolics.gradient(L, z)

function residual(w, θ, κ)
    # primals 
    z = w[1:10]
    pa, pb, ta, tb, sa, sb = unpack_primals(z)

    # duals 
    y = w[10 .+ (1:12)]
    λa, λb, νa, νb = unpack_duals(y)

    # Lagrangian 
    lag = lagrangian(z, y, θ) #+ 1.0e-2 * dot(λa, λa) - 1.0e-2 * dot(νa, νa) - 1.0e-5 * dot(νb, νb)
    lagz = Symbolics.gradient(lag, z)
    con = constraints(z, θ) 
    comp = [
            ta[1] * νa[1]; 
            sa[1] * νa[2]; 
            tb[1] * νb[1];
            sb[1] * νb[2];
           ]

    res = [
            lagz;
            con;
            comp .- κ;
          ]

    return res 
end

w = [z; y]
r = residual([z; y], θ, κ)
rw = Symbolics.jacobian(r, [z; y])
rθ = Symbolics.jacobian(r, θ)

r_func = eval(Symbolics.build_function(r, w, θ, κ)[2])
rw_func = eval(Symbolics.build_function(rw, w, θ)[2])
rθ_func = eval(Symbolics.build_function(rθ, w, θ)[2])

# pre-allocate
r0 = zeros(nw)
rw0 = zeros(nw, nw)
rθ0 = zeros(nw, nθ)

# random 
# w0 = randn(nw) 
# θ0 = randn(nθ) 
# κ0 = randn(1)

# r_func(r0, w0, θ0, κ0)
# rw_func(rw0, w0, θ0)
# rθ_func(rθ0, w0, θ0)

function rw_func_reg(rw, w, θ) 
    rw_func(rw, w, θ)
    rw .+= Diagonal([1.0e-5 * ones(10); -1.0e-5 * ones(8); -1.0e-5 * ones(4)])
    return 
end


## setup 
xa = [-1.0; 0.0; 0.0]
qa = Quaternion(RotY(0.5 * π) * RotX(0.0 * π))
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
w0 = [xa; xb; 0.5; 0.5; 1.0; 1.0; zeros(8); ones(4)]
θ0 = [pa1; pa2; pb1; pb2]
κ0 = [1.0]

r_func(r0, w0, θ0, κ0)
rw_func_reg(rw0, w0, θ0)
# rw_func(rw0, w0, θ0)
rθ_func(rθ0, w0, θ0)
cond(rw0)
rank(rw0)
rw0 \ r0
rw0 \ rθ0

## solver 
idx = RoboDojo.IndicesOptimization(
    nw, nw,
    [collect(7:10), collect(19:22)], [collect(7:10), collect(19:22)],
    Vector{Vector{Vector{Int}}}(), Vector{Vector{Vector{Int}}}(),
    collect(1:18),
    collect(19:22),
    Vector{Int}(),
    Vector{Vector{Int}}(),
    collect(19:22),
)

ip = RoboDojo.interior_point(w0, θ0;
    s = RoboDojo.Euclidean(length(w0)),
    idx = idx,
    r! = r_func, rz! = rw_func_reg, rθ! = rθ_func,
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


pa = ip.z[0 .+ (1:3)]
pb = ip.z[3 .+ (1:3)]
dir = (pa - pb) ./ norm(pa - pb)

ip.z[6 .+ (1:1)]
ip.z[7 .+ (1:1)]
ip.z[8 .+ (1:1)]
ip.z[9 .+ (1:1)]

ip.z[18 .+ (1:4)]

# ∂pa∂xa = ip.δz[0 .+ (1:3), 0  .+ (1:3)]
# ∂pa∂qa = ip.δz[0 .+ (1:3), 3  .+ (1:4)]
# ∂pa∂xb = ip.δz[0 .+ (1:3), 7  .+ (1:3)]
# ∂pa∂qb = ip.δz[0 .+ (1:3), 10 .+ (1:4)]

# ∂pb∂xa = ip.δz[3 .+ (1:3), 0  .+ (1:3)]
# ∂pb∂qa = ip.δz[3 .+ (1:3), 3  .+ (1:4)]
# ∂pb∂xb = ip.δz[3 .+ (1:3), 7  .+ (1:3)]
# ∂pb∂qb = ip.δz[3 .+ (1:3), 10 .+ (1:4)]

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
settransform!(vis[:cap1], compose(Translation(xa), LinearMap(Quaternion(qa...)))
settransform!(vis[:cap2], compose(Translation(xb), LinearMap(Quaternion(qb...)))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(pa - dir * ra))
settransform!(vis[:cs2], Translation(pb + dir * rb))