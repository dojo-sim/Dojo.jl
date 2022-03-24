# https://wickedengine.net/2020/04/26/capsule-collision-detection/

struct CapsuleContact{T}
    position::Vector{T} 
    orientation::Quaternion{T}
    height::T 
    radius::T 
end

p1 = [0.0, 0.0, 0.0] 
# o1 = one(Quaternion)
o1 = Quaternion(RotX(0.6) * RotY(0.0 * π))
h1 = 1.0 
r1 = 0.1 
cap1 = CapsuleContact(p1, o1, h1, r1)

p2 = [1.0, 0.0, 0.0] 
o2 = one(Quaternion)
o2 = Quaternion(RotZ(0.0) * RotY(-0.25 * π))
h2 = 1.0 
r2 = 0.1 
cap2 = CapsuleContact(p2, o2, h2, r2)

function get_tip_base(cap::CapsuleContact) 
    k = [0.0, 0.0, 0.5 * cap.height + cap.radius]
    tip = cap.position + vector_rotate(k, cap.orientation)
    base = cap.position + vector_rotate(-k, cap.orientation)
    return Array(tip), Array(base) 
end

t1, b1 = get_tip_base(cap1)
t2, b2 = get_tip_base(cap2)

function closest_point_on_line(pa::Vector{T}, pb::Vector{T}, p::Vector{T}) where T 
    pab = pb - pa 
    t = dot(p - pa, pab) / dot(pab, pab)
    return pa + min(max(t, 0.0), 1.0) * pab
end

function closest_segment_points_capsule_capsule(capa::CapsuleContact, capb::CapsuleContact)
    # capsule a
    a_tip, a_base = get_tip_base(capa)
    a_radius = capa.radius

    # capsule b
    b_tip, b_base = get_tip_base(capb)
    b_radius = capb.radius

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
    # # direction 
    # direction = bestB - bestA
    # direction_magnitude = norm(direction)
    # direction ./= direction_magnitude

    # c1 = bestA + direction * a_radius
    # c2 = bestB - direction * b_radius
    # # distance = direction_magnitude - (a_radius + b_radius)
    # return c1, c2
end


function kinematics_capsule_capsule(capa::CapsuleContact, capb::CapsuleContact)
    bestA, bestB = closest_segment_points_capsule_capsule(capa, capb)

    kA = vector_rotate(bestA - capa.position, inv(capa.orientation))
    kB = vector_rotate(bestB - capb.position, inv(capb.orientation))

    return kA, kB
end

ka, kb = kinematics_capsule_capsule(cap1, cap2)

function closest_points_capsule_capsule(capa::CapsuleContact, capb::CapsuleContact)
    @show "hi"
    bestA, bestB = closest_segment_points_capsule_capsule(capa, capb)

    direction = bestB - bestA
    direction_magnitude = norm(direction)
    direction ./= direction_magnitude

    c1 = bestA + direction * capa.radius
    c2 = bestB - direction * capb.radius
    # distance = direction_magnitude - (a_radius + b_radius)
    return c1, c2
end

function direction_capsule_capsule(capa::CapsuleContact, capb::CapsuleContact)
    ca, cb = closest_segment_points_capsule_capsule(capa, capb)

    direction = cb - ca
    direction_magnitude = norm(direction)

    if direction_magnitude > 0.0
        direction ./= direction_magnitude
    else
        direction = [1.0, 0.0, 0.0]
    end

    return direction 
end

function sdf_capsule_capsule(capa::CapsuleContact, capb::CapsuleContact)
    ca, cb = closest_points_capsule_capsule(capa, capb)

    direction = cb - ca
    direction_magnitude = norm(direction)

    return direction_magnitude
end


c1, c2 = closest_points_capsule_capsule(cap1, cap2)

# dir = direction_capsule_capsule(cap1, cap2)

# dist = distance_capsule_capsule(cap1, cap2)

# pa = [0.0, 0.0, 0.0]
# pb = [1.0, 1.0, 0.0]
# p = [1.0, 0.5, 0.0]

# closest_point_on_line(pa, pb, p)

# # capsule a
# a_tip = [1.0, 1.0, 0.0]
# a_base = [0.0, 0.0, 0.0]
# a_radius = 0.01

# # capsule b
# b_tip = [2.0, 1.0, 0.0]
# b_base = [2.0, 0.0, 0.0]
# b_radius = 0.01

# # capsule a 
# a_normal = a_tip - a_base 
# a_normal ./= norm(a_normal)
# a_line_end_offset = a_normal * a_radius
# a_A = a_base + a_line_end_offset 
# a_B = a_tip - a_line_end_offset 

# # capsule b 
# b_normal = b_tip - b_base 
# b_normal ./= norm(b_normal) 
# b_line_end_offset = b_normal * b_radius 
# b_A = b_base + b_line_end_offset 
# b_B = b_tip - b_line_end_offset

# # vectors between line endpoints 
# v0 = b_A - a_A 
# v1 = b_B - a_A 
# v2 = b_A - a_B 
# v3 = b_B - a_B

# # squared distances 
# d0 = dot(v0, v0) 
# d1 = dot(v1, v1) 
# d2 = dot(v2, v2) 
# d3 = dot(v3, v3) 

# # best potential endpoint on capsule a
# if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1)
#     bestA = a_B
# else
#     bestA = a_A
# end
# # select point on capsule B line segment nearest best potential endpoint on A capsule
# bestB = closest_point_on_line(b_A, b_B, bestA)

# # same for capsule A segment
# bestA = closest_point_on_line(a_A, a_B, bestB)

# # direction 
# direction = bestB - bestA
# direction_magnitude = norm(direction)
# direction ./= direction_magnitude

# distance = direction_magnitude - (a_radius + b_radius)

# visualize 
vis = Visualizer()
open(vis)

# capsule a
cyl1 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * cap1.height), Point(0.0, 0.0, 0.5 * cap1.height), cap1.radius)
sph1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), cap1.radius)
color1 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap1][:cyl], cyl1, MeshPhongMaterial(color=color1))
setobject!(vis[:cap1][:base], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:base], Translation(0.0, 0.0, -0.5 * cap1.height))
setobject!(vis[:cap1][:tip], sph1, MeshPhongMaterial(color=color1))
settransform!(vis[:cap1][:tip], Translation(0.0, 0.0, 0.5 * cap1.height))

# capsule b
cyl2 = GeometryBasics.Cylinder(Point(0.0, 0.0, -0.5 * cap2.height), Point(0.0, 0.0, 0.5 * cap2.height), cap2.radius)
sph2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), cap2.radius)
color2 = Colors.RGBA(0.7, 0.7, 0.7, 0.5);
setobject!(vis[:cap2][:cyl], cyl2, MeshPhongMaterial(color=color2))
setobject!(vis[:cap2][:base], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:base], Translation(0.0, 0.0, -0.5 * cap2.height))
setobject!(vis[:cap2][:tip], sph2, MeshPhongMaterial(color=color2))
settransform!(vis[:cap2][:tip], Translation(0.0, 0.0, 0.5 * cap2.height))

# set configuration
settransform!(vis[:cap1], compose(Translation(cap1.position), LinearMap(cap1.orientation)))
settransform!(vis[:cap2], compose(Translation(cap2.position), LinearMap(cap2.orientation)))

# closest points
color_cp = Colors.RGBA(0.0, 0.0, 0.0, 1.0);
cs1 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
cs2 = GeometryBasics.Sphere(Point(0.0, 0.0, 0.0), 0.025)
setobject!(vis[:cs1], cs1, MeshPhongMaterial(color=color_cp))
setobject!(vis[:cs2], cs2, MeshPhongMaterial(color=color_cp))
settransform!(vis[:cs1], Translation(c1))
settransform!(vis[:cs2], Translation(c2))

