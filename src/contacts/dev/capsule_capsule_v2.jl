# https://wickedengine.net/2020/04/26/capsule-collision-detection/

p1 = [0.0, 0.0, 0.0] 
o1 = one(UnitQuaternion)
h1 = 1.0 
r1 = 0.1 

p2 = [1.0, 0.0, 0.0] 
o2 = one(UnitQuaternion)
# o2 = UnitQuaternion(RotZ(0.41) * RotY(-0.3 * π))
h2 = 1.0 
r2 = 0.1 

function get_tip_base(position, orientation, height, radius) 
    k = [0.0, 0.0, 0.5 * height + radius]
    tip = position + Array(vector_rotate(k, orientation))
    base = position + Array(vector_rotate(-k, orientation))
    return tip, base
end

t1, b1 = get_tip_base(p1, o1, h1, r1)
t2, b2 = get_tip_base(p2, o2, h2, r2)

function closest_point_on_line(pa::Vector{T}, pb::Vector{T}, p::Vector{T}) where T 
    pab = pb - pa 
    t = dot(p - pa, pab) / dot(pab, pab)
    return pa + min(max(t, 0.0), 1.0) * pab
end

function closest_segment_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    # capsule a
    a_tip, a_base = get_tip_base(xa, qa, ha, ra)
    a_radius = ra

    # capsule b
    b_tip, b_base = get_tip_base(xb, qb, hb, rb)
    b_radius = rb

    # capsule a 
    a_normal = a_tip - a_base 
    a_normal ./= norm(a_normal)
    a_line_end_offset = a_normal * a_radius
    a_A = a_base + a_line_end_offset 
    a_B = a_tip - a_line_end_offset 

    # capsule b 
    b_normal = b_tip - b_base 
    b_normal ./= norm(b_normal) 
    b_line_end_offset = b_normal * b_radius 
    b_A = b_base + b_line_end_offset 
    b_B = b_tip - b_line_end_offset

    # vectors between line endpoints 
    v0 = b_A - a_A 
    v1 = b_B - a_A 
    v2 = b_A - a_B 
    v3 = b_B - a_B

    # squared distances 
    d0 = dot(v0, v0) 
    d1 = dot(v1, v1) 
    d2 = dot(v2, v2) 
    d3 = dot(v3, v3) 

    # best potential endpoint on capsule a
    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1)
        bestA = a_B
    else
        bestA = a_A
    end

    # select point on capsule B line segment nearest best potential endpoint on A capsule
    bestB = closest_point_on_line(b_A, b_B, bestA)

    # same for capsule A segment
    bestA = closest_point_on_line(a_A, a_B, bestB)

    return bestA, bestB
end


function kinematics_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    bestA, bestB = closest_segment_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)

    kA = vector_rotate(bestA - xa, inv(qa))
    kB = vector_rotate(bestB - xb, inv(qb))

    return kA, kB
end

ka, kb = kinematics_capsule_capsule(p1, o1, h1, r1, p2, o2, h2, r2)

function closest_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    bestA, bestB = closest_segment_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)

    direction = bestB - bestA
    direction_magnitude = norm(direction)
    direction ./= direction_magnitude

    c1 = bestA + direction * ra
    c2 = bestB - direction * rb
    return c1, c2
end

function direction_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    ca, cb = closest_segment_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)

    direction = cb - ca
    direction_magnitude = norm(direction)

    if direction_magnitude > 0.0
        direction ./= direction_magnitude
    else
        direction = [1.0, 0.0, 0.0]
    end

    return direction 
end

function sdf_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)
    ca, cb = closest_points_capsule_capsule(xa, qa, ha, ra, xb, qb, hb, rb)

    direction = cb - ca
    direction_magnitude = norm(direction)

    return direction_magnitude
end

c1, c2 = closest_points_capsule_capsule(p1, o1, h1, r1, p2, o2, h2, r2)
sdf_capsule_capsule(p1, o1, h1, r1, p2, o2, h2, r2)

# visualize 
vis = Visualizer()
open(vis)

# capsule a
cyl1 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * h1), Point(0.0, 0.0, 0.5 * h1), r1)
sph1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), r1)
color1 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap1][:cyl], cyl1, MeshPhongMaterial(color=color1))
setobject!(vis[:cap1][:base], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:base], Translation(0.0, 0.0, -0.5 * h1))
setobject!(vis[:cap1][:tip], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:tip], Translation(0.0, 0.0, 0.5 * h1))

# capsule b
cyl2 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * h2), Point(0.0, 0.0, 0.5 * h2), r2)
sph2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), r2)
color2 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap2][:cyl], cyl2, MeshPhongMaterial(color=color2))
setobject!(vis[:cap2][:base], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:base], Translation(0.0, 0.0, -0.5 * h2))
setobject!(vis[:cap2][:tip], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:tip], Translation(0.0, 0.0, 0.5 * h2))

# set configuration
settransform!(vis[:cap1], compose(Translation(p1), LinearMap(o1)))
settransform!(vis[:cap2], compose(Translation(p2), LinearMap(o2)))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(c1))
settransform!(vis[:cs2], Translation(c2))

