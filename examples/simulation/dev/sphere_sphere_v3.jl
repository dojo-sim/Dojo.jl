using LinearAlgebra
using SparseArrays

using MeshCat
using GeometryBasics
using Colors
using CoordinateTransformations
using FiniteDiff
using Rotations

# visualize
function set_background!(vis::Visualizer; top_color=RGBA(1,1,1.0), bottom_color=RGBA(1,1,1.0))
    setprop!(vis["/Background"], "top_color", top_color)
    setprop!(vis["/Background"], "bottom_color", bottom_color)
end
vis = Visualizer()
open(vis)
# set_background!(vis)

## setup 
xa = [0.0; 0.0; 0.0]
qa = [1.0; 0.0; 0.0; 0.0]
ra = 0.25

xb = [1.0; 1.0e-13; 0.0]
qb = [1.0; 0.0; 0.0; 0.0]
rb = 0.25
collision = SphereSphereCollision{Float64,3,3,9}(
        szeros(3),
        szeros(3),
        ra, 
        rb)

pa = contact_point(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))
pb = contact_point(:child,  collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))

n = contact_normal(collision,  xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))
t = contact_tangent(collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))

# sphere a
sa = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), ra)
color1 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:spherea], sa, MeshPhongMaterial(color=color1))
settransform!(vis[:spherea], Translation(xa...))

# sphere b
sb = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), rb)
color2 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:sphereb], sb, MeshPhongMaterial(color=color2))
settransform!(vis[:sphereb], Translation(xb...))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(pa...))
settransform!(vis[:cs2], Translation(pb...))

# normals (NOTE: visualize opposite directions)
set_arrow!(vis, pa, -vec(n), name=:parent_normal)
set_arrow!(vis, pb, vec(n), name=:child_normal)

set_arrow!(vis, pa, vec(t[1, :]), name=:parent_tangent1, color=RGBA(0.0, 1.0, 0.0, 1.0))
set_arrow!(vis, pa, vec(t[2, :]), name=:parent_tangent2, color=RGBA(0.0, 0.0, 1.0, 1.0))

set_arrow!(vis, pb, -1.0 * vec(t[1, :]), name=:child_tangent1, color=RGBA(0.0, 1.0, 0.0, 1.0))
set_arrow!(vis, pb, -1.0 * vec(t[2, :]), name=:child_tangent2, color=RGBA(0.0, 0.0, 1.0, 1.0))




## Jacobians contact_normal(collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))


# FiniteDiff.finite_difference_jacobian(x -> vec(contact_normal(collision, x, Quaternion(qa..., false), xb, Quaternion(qb..., false))), xa)

λ = rand(3)
λ' * FiniteDiff.finite_difference_jacobian(x -> contact_normal(collision, x, Quaternion(qa..., false), xb, Quaternion(qb..., false)), xa)
λ' * FiniteDiff.finite_difference_jacobian(q -> contact_normal(collision, xa, Quaternion(q..., false), xb, Quaternion(qb..., false)), vector(qa))

∂contact_normal_jvp∂x(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
∂contact_normal_jvp∂q(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)

contact_tangent(collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false))

λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, Quaternion(qa..., false), xb, Quaternion(qb..., false))[1, :]', xa)
λ' * FiniteDiff.finite_difference_jacobian(x -> contact_tangent(collision, x, Quaternion(qa..., false), xb, Quaternion(qb..., false))[2, :]', xa)

λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xa, Quaternion(q..., false), xb, Quaternion(qb..., false))[1, :], vector(qa))
λ' * FiniteDiff.finite_difference_jacobian(q -> contact_tangent(collision, xa, Quaternion(q..., false), xb, Quaternion(qb..., false))[2, :], vector(qa))

∂contact_tangent_jvp∂x(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
∂contact_tangent_jvp∂q(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)

∂contact_tangent_jvp∂x(:child, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
∂contact_tangent_jvp∂q(:child, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)

λ = rand(2)
∂contact_tangent_vjp∂x(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
∂contact_tangent_vjp∂q(:parent, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)

∂contact_tangent_vjp∂x(:child, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
∂contact_tangent_vjp∂q(:child, collision, xa, Quaternion(qa..., false), xb, Quaternion(qb..., false), λ)
