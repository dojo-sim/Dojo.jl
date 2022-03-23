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

## setup 
xa = [-1.0; 0.0; 0.0]
qa = Quaternion(RotY(0.0 * π) * RotX(0.0 * π))
ra = 0.1
ha = 0.2

xb = [1.0; 0.0; 0.0]
qb = Quaternion(RotZ(0.0 * π) * RotY(0.0 * π) * RotX(0.0))
rb = 0.1 
hb = 0.2

collision = CapsuleCapsuleCollision(ra, ha, rb, hb)

pa = capsule_contact_origin(:parent, collision, xa, qa, xb, qb)
pb = capsule_contact_origin(:child, collision, xa, qa, xb, qb)
n = contact_normal(collision, xa, qa, xb, qb)
dir = vec(n)

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
settransform!(vis[:cap1], compose(Translation(xa), LinearMap(Quaternion(qa))))
settransform!(vis[:cap2], compose(Translation(xb), LinearMap(Quaternion(qb))))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(pa - dir * ra))
settransform!(vis[:cs2], Translation(pb + dir * rb))