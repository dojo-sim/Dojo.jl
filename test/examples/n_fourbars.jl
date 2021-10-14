using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x, y = .1, .1

vert11 = [0.;0.;l1 / 2]
vert12 = -vert11
verts = [[vert11];[vert12]]

# Initial orientation
offset1 = pi / 4
offset2 = pi / 2
phi1 = pi / 8
q1 = UnitQuaternion(RotX(phi1))
qoff1 = UnitQuaternion(RotX(offset1))
qoff2 = UnitQuaternion(RotX(offset2))

N = 10

# Links
origin = Origin{Float64}()
links = [Box(x, y, l1, l1, color = RGBA(1., 1., 0.)) for i = 1:4 * N]


function fourbar(links, vertices, axis)
    j1 = EqualityConstraint(Revolute(links[1], links[2], axis; p1=vertices[1], p2=vertices[2]))
    j2 = EqualityConstraint(Revolute(links[2], links[3], axis; p1=vertices[3], p2=vertices[2]), Cylindrical(links[2], links[4], axis; p1=vertices[2], p2=vertices[2]))
    j3 = EqualityConstraint(Revolute(links[4], links[5], axis; p1=vertices[3], p2=vertices[2]))
    j4 = EqualityConstraint(Revolute(links[3], links[5], axis; p1=vertices[3], p2=vertices[3]))

    return j1, j2, j3, j4
end

function initfourbar!(mechanism, links, vertices, Δq1, Δq2)
    setPosition!(links[1], links[2], p1 = vertices[1], p2 = vertices[2], Δq = Δq1)
    setPosition!(links[2], links[3], p1 = vertices[3], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(links[2], links[4], p1 = vertices[2], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(links[4], links[5], p1 = vertices[3], p2 = vertices[2], Δq = Δq2 * Δq2)

    return 
end

# Constraints
constraints = [fourbar([origin;links[1:4]], [[zeros(3)];verts], ex)...]
for i = 2:N
    push!(constraints, fourbar(links[(i - 1) * 4:i * 4], [[vert12];verts], ex)...)
end


mech = Mechanism(origin, links, constraints)
if N > 1
    initfourbar!(mech, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff1, q1)
else
    initfourbar!(mech, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff2, q1)
end

for i = 2:N - 1
    initfourbar!(mech, links[(i - 1) * 4:i * 4], [[vert12];verts], q1 / q1, q1)
end

if N > 1
    initfourbar!(mech, links[(N - 1) * 4:N * 4], [[vert12];verts], inv(qoff1) * qoff2, q1)
end

