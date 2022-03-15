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
"""
    variables:
    pa ∈ R^3 # parent body contact point
    pb ∈ R^3 # child  body  contact point 

    sa ∈ R+  # slack for parent body collision set 
    sb ∈ R+  # slack for child  body collision set

    λa ∈ R   # dual  for parent equality constraint 
    λb ∈ R   # dual  for child  equality constraint 

    νa ∈ R+  # dual  for parent slack 
    νb ∈ R+  # dual  for child  slack 
    
    data: 
    xa ∈ R^3 # parent body position 
    qa ∈ H   # parent body orientation 
    xb ∈ R^3 # child  body position 
    qb ∈ H   # child  body orientation 
    ra ∈ R+  # parent body collision radius
    rb ∈ R+  # child  body collision radius

    objective: 
    ||pa - pb||

    -> (pa - pb)' (pa - pb) # squared distance

    constraints: 
    ||R(qa)' * (pa - xa)||_Σa^-1 <= 1.0 # ellipsoid
    ||R(qb)' * (pb - xb)||_Σb^-1 <= 1.0 # ellipsoid
        where Σ^-1 = Diagonal([1.0 / x^2, 1.0 / y^2, 1.0 / z^2])
    -> 
        sa - (1.0 - (pa - xa)' R(qa) Σa^-1 R(qa)' (pa - xa)) = 0 # squared + slack
        sb - (1.0 - (pb - xb)' R(qb) Σb^-1 R(qb)' (pb - xb)) = 0 # squared + slack
        sa, sb >= 0

    primal: 
    z = (pa, pb, sa, sb) ∈ R^8
    dual:
    y = (λa, λb, νa, νb) ∈ R^4
    data:
    θ = (xa, qa, xb, qb, ax, ay, az, bx, by, bz) ∈ R^16
"""

nz = 12
ny = 12
nw = nz + ny 
nθ = 18
@variables z[1:nz] y[1:ny] θ[1:nθ] κ[1:1]

function unpack_primals(z) 
    pa = z[0 .+ (1:3)]
    pb = z[3 .+ (1:3)]
    sa = z[6 .+ (1:3)]
    sb = z[9 .+ (1:3)]
    return pa, pb, sa, sb 
end

function unpack_duals(y)
   λa = y[0 .+ (1:3)]
   λb = y[3 .+ (1:3)]
   νa = y[6 .+ (1:3)]
   νb = y[9 .+ (1:3)]

   return λa, λb, νa, νb
end

function unpack_data(θ)
    xa = θ[1:3]
    qa = θ[3 .+ (1:4)]
    xb = θ[7 .+ (1:3)]
    qb = θ[10 .+ (1:4)]
    ra = θ[14 .+ (1:1)]
    ha = θ[15 .+ (1:1)]
    rb = θ[16 .+ (1:1)]
    hb = θ[17 .+ (1:1)]
    
    return xa, qa, xb, qb, ra, ha, rb, hb
end

function objective(z)
    pa, pb, sa, sb = unpack_primals(z)
    d = pa - pb 
    return dot(d, d)
end

function constraints(z, θ) 
    pa, pb, sa, sb = unpack_primals(z)
    xa, qa, xb, qb, ra, ha, rb, hb = unpack_data(θ)

    da = pa - xa 
    db = pb - xb 

    [
        sa[1] - (ra[1] - norm(da[1:2]));
        sa[2] - (0.5 * ha[1] + sa[1] - da[3]);
        sa[3] - (da[3] + 0.5 * ha[1] + sa[1]);

        sb[1] - (rb[1] - norm(db[1:2]));
        sb[2] - (0.5 * hb[1] + sb[1] - db[3]);
        sb[3] - (db[3] + 0.5 * hb[1] + sb[1]);
    ]
end

function lagrangian(z, y, θ)
    # initialize
    L = 0.0

    # primals 
    pa, pb, sa, sb = unpack_primals(z)

    # duals 
    λa, λb, νa, νb = unpack_duals(y)

    # objective 
    J = objective(z) 
    L += J 

    # constraints
    c = constraints(z, θ)
    L += sum([λa[1] λa[2] λa[3] λb[1] λb[2] λb[3]] * c)

    # inequalities 
    L -= sa[1] * νa[1]
    L -= sa[2] * νa[2]
    L -= sa[3] * νa[3]

    L -= sb[1] * νb[1] 
    L -= sb[2] * νb[2] 
    L -= sb[3] * νb[3] 

    return L 
end

L = lagrangian(z, y, θ)
Lz = Symbolics.gradient(L, z)

function residual(w, θ, κ)
    # primals 
    z = w[1:12]
    pa, pb, sa, sb = unpack_primals(z)

    # duals 
    y = w[12 .+ (1:12)]
    λa, λb, νa, νb = unpack_duals(y)

    # Lagrangian 
    lag = lagrangian(z, y, θ)
    lagz = Symbolics.gradient(lag, z)
    con = constraints(z, θ) 
    comp = [
            sa[1] * νa[1]; 
            sa[2] * νa[2]; 
            sa[3] * νa[3]; 
            sb[1] * νb[1];
            sb[2] * νb[2];
            sb[3] * νb[3];
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
    rw .+= Diagonal([1.0e-5 * ones(12); -1.0e-5 * ones(12)])
    return 
end

# (rw0 + Diagonal([1.0e-5 * ones(8); -1.0e-5 * ones(4)])) \ rθ0

## setup 
xa = [-1.0; 0.0; 0.0]
qa = UnitQuaternion(RotY(0.0) * RotX(0.0 * π))
qa = [qa.w; qa.x; qa.y; qa.z]
ra = 0.1
ha = 0.2

xb = [1.0; 0.0; 0.0]
qb = UnitQuaternion(RotZ(0.0 * π) * RotY(0.0) * RotX(0.0))
qb = [qb.w; qb.x; qb.y; qb.z]
rb = 0.1 
hb = 0.2

## initialization 
w0 = [xa + 1.0e-1 * randn(3); xb + 1.0e-1 * randn(3); ones(6)...; zeros(6)...; ones(6)...]
θ0 = [xa; qa; xb; qb; ra; ha; rb; hb]
κ0 = [1.0]

r_func(r0, w0, θ0, κ0)
rw_func_reg(rw0, w0, θ0)
# rw_func(rw0, w0, θ0)
rθ_func(rθ0, w0, θ0)
cond(rw0)
rw0 \ rθ0

## solver 
idx = RoboDojo.IndicesOptimization(
    nw, nw,
    [collect(7:12), collect(19:24)], [collect(7:12), collect(19:24)],
    Vector{Vector{Vector{Int}}}(), Vector{Vector{Vector{Int}}}(),
    collect(1:18),
    collect(19:24),
    Vector{Int}(),
    Vector{Vector{Int}}(),
    collect(19:24),
)

ip = RoboDojo.interior_point(w0, θ0;
    s = RoboDojo.Euclidean(length(w0)),
    idx = idx,
    r! = r_func, rz! = rw_func_reg, rθ! = rθ_func,
    r  = zeros(idx.nΔ),
    rz = zeros(idx.nΔ, idx.nΔ),
    rθ = zeros(idx.nΔ, length(θ0)),
    opts = RoboDojo.InteriorPointOptions(
            undercut=10.0,
            γ_reg=0.1,
            r_tol=1e-6,
            κ_tol=1e-6,  
            max_ls=50,
            ϵ_min=0.01,
            diff_sol=true,
            verbose=true))

RoboDojo.interior_point_solve!(ip)

pa = ip.z[0 .+ (1:3)]
pb = ip.z[3 .+ (1:3)]

∂pa∂xa = ip.δz[0 .+ (1:3), 0  .+ (1:3)]
∂pa∂qa = ip.δz[0 .+ (1:3), 3  .+ (1:4)]
∂pa∂xb = ip.δz[0 .+ (1:3), 7  .+ (1:3)]
∂pa∂qb = ip.δz[0 .+ (1:3), 10 .+ (1:4)]

∂pb∂xa = ip.δz[3 .+ (1:3), 0  .+ (1:3)]
∂pb∂qa = ip.δz[3 .+ (1:3), 3  .+ (1:4)]
∂pb∂xb = ip.δz[3 .+ (1:3), 7  .+ (1:3)]
∂pb∂qb = ip.δz[3 .+ (1:3), 10 .+ (1:4)]

# sphere a
sa = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 1.0)
color1 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:spherea], sa, MeshPhongMaterial(color=color1))
settransform!(vis[:spherea], compose(Translation(xa...), LinearMap(UnitQuaternion(qa...) * Diagonal([ax, ay, az]))))

# sphere b
sb = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 1.0)
color2 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:sb], sb, MeshPhongMaterial(color=color2))
settransform!(vis[:sb], compose(Translation(xb...), LinearMap(UnitQuaternion(qb...) * Diagonal([bx, by, bz]))))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(pa...))
settransform!(vis[:cs2], Translation(pb...))


