# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))

# Open visualizer
vis = Visualizer()
open(vis)

# Build mechanism
include("mechanism_zoo.jl")


###############################################################################
# METHODS NEW EQUALITY CONSTRAINT
################################################################################

# t2r3
function ForcePrismatic(body1::AbstractBody{T}, body2, axis;
        p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}),
        spring = zero(T), damper = zero(T)) where T
    Translational2{T}(body1, body2; p1, p2, axis, spring, damper),
    Force121{T}(body1, body2; p1, p2, axis, spring, damper),
    Rotational3{T}(body1, body2; qoffset, spring, damper)
end

################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################

# Parameters
ex = [0; 0; 1.0]
h = 1.
r = .05
vert11 = [0; r; 0.0]
vert12 = -vert11
Nlink = 2

# Links
origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

# Constraints
jointb1 = EqualityConstraint(Fixed(origin, links[1]; p1 = zeros(3), p2 = zeros(3)))
if Nlink > 1
    eqcs = [
        jointb1;
        [EqualityConstraint(ForcePrismatic(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = 1.0, damper = 0.1)) for i = 2:Nlink]
        ]
else
    eqcs = [jointb1]
end
mech = Mechanism(origin, links, eqcs, g = -9.81, Δt = 0.01)

# mech.eqconstraints[1]
# eqc2 = mech.eqconstraints[2]
# eqc2.λsol
# tra2 = eqc2.constraints[1]
# force2 = eqc2.constraints[2]
# rot2 = eqc2.constraints[3]
# springforce(force2, body1.state, body2.state)
# springforce()

# mech = getmechanism(:nslider, Nlink = 5)
initialize!(mech, :nslider)
# storage = simulate!(mech, 1.0, record = true, solver = :mehrotra!)
storage = simulate!(mech, 10.0, record = true, solver = :mehrotra!)

visualize(mech, storage, vis = vis)




eqc1 = mech.eqconstraints[1]
eqc2 = mech.eqconstraints[2]
tra1 = eqc1.constraints[1]
force1 = eqc1.constraints[2]
rot1 = eqc1.constraints[3]
force2 = eqc2.constraints[2]
body1, body2 = mech.bodies
origin
body1
body2

eqc1.λinds
∂g∂ʳposa(tra1, body1, body2, body2.id)


constraintmat(tra1)
A = constraintmat(force1)
Aᵀ = zerodimstaticadjoint(A)
Aᵀ * A * ([0.1] * A)

([0.1] * A)
∂g∂ʳvelb(mech, eqc1, body1)
∂g∂ʳvelb(tra1, body1, body2, body2.id, mech.Δt)
∂g∂ʳvelb(force1, body1, body2, body2.id, mech.Δt)
∂g∂ʳvelb(rot1, body1, body2, body2.id, mech.Δt)




∂g∂ʳself(mech, eqc1)
∂g∂ʳself(tra1)
∂g∂ʳself(force1)
∂g∂ʳself(rot1)
g(mech, eqc1)
g(mech, eqc2)
gc(mech, eqc2)

constraintmat(force1)
g(tra1, origin, body1, mech.Δt, rand(33))
g(rot1, origin, body1, mech.Δt, rand(33))
g(force1, origin, body1)
g(force1, body1, body2)


λindsf1 = Vector(eqc1.λinds[2][1]:eqc1.λinds[2][2])
eqc1.λsol[2][λindsf1]
g(force1, origin, body1)
g(force1, body1, body2)

springforce(force1, body1.state)
springforce(force1, body1.state, body2.state)

damperforce(force1, body1.state)
damperforce(force1, body1.state, body2.state)


eqc1.indss
eqc2.inds

Vector(SUnitRange(eqc1.inds[2][1], eqc1.inds[2][2]))



tra, ftr, rot = ForcePrismatic(origin, links[1], ex; p2 = vert11, spring = spring0, damper = damper0)
tra = tra[1]
ftr = ftr[1]
rot = rot[1]

constraintmat(tra)
constraintmat(ftr)
constraintmat(rot)
nullspacemat(tra)
nullspacemat(ftr)
nullspacemat(rot)




jointb1.inds
eqcs[1].inds
eqcs[2].inds
