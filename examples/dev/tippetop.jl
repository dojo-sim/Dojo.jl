# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))



function gettippetop(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, contact::Bool=true, contact_mode::Symbol = :soc) where T
    origin = Origin{T}(name="origin")
    radius = 0.5
    mass = 1.0
    α = 0.2
    bodies = [Sphere(radius, mass, name="sphere1"), Sphere(radius*α, mass*α^3, name="sphere2")]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = "floating_joint"),
        EqualityConstraint(Fixed(bodies[1], bodies[2], p1=[0,0,radius], p2=zeros(3)), name = "fixed_joint"),]
    mechanism = Mechanism(origin, bodies, eqcs, Δt = Δt, g = g)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        (contact_mode == :soc) && (contineqcs = [
            contactconstraint(getbody(mechanism, "sphere1"), normal, cf, p=contact, offset=[0,0,radius]),
            contactconstraint(getbody(mechanism, "sphere2"), normal, cf, p=contact, offset=[0,0,radius*α])
            ])
        (contact_mode == :linear) && (contineqcs = [
            linearcontactconstraint(getbody(mechanism, "sphere1"), normal, cf, p=contact, offset=[0,0,radius]),
            linearcontactconstraint(getbody(mechanism, "sphere2"), normal, cf, p=contact, offset=[0,0,radius*α])
            ])
        (contact_mode == :impact) && (contineqcs = [
            impactconstraint(getbody(mechanism, "sphere1"), normal, p=contact, offset=[0,0,radius]),
            impactconstraint(getbody(mechanism, "sphere2"), normal, p=contact, offset=[0,0,radius*α])
            ])
        setPosition!(mechanism, geteqconstraint(mechanism, "floating_joint"), [0;0;radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, eqcs, contineqcs, g=g, Δt=Δt)
    end
    return mechanism
end

function initializetippetop!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}

    eqc2 = geteqconstraint(mechanism, "fixed_joint")
    radius = eqc2.constraints[1].vertices[1][3]
    origin = mechanism.origin
    body1 = getbody(mech, "sphere1")
    body2 = getbody(mech, "sphere2")

    zeroVelocity!(mechanism)
    # setPosition!(mechanism, eqc, [x; rotation_vector(q)])
    # setVelocity!(mechanism, eqc, [v; ω])
    setPosition!(origin, body1; p1 = [0;0;radius], p2 = [0;0;0], Δx = x, Δq = q)
    setPosition!(body1,  body2; p1 = [0;0;radius], p2 = [0;0;0], Δx = [0;0;0], Δq = one(UnitQuaternion))
    setVelocity!(origin, body1; p1 = [0;0;radius], p2 = [0;0;0], Δv = v, Δω = ω)
    setVelocity!(body1,  body2; p1 = [0;0;radius], p2 = [0;0;0], Δv = [0;0;0], Δω = [0;0;0])
    return nothing
end





mech = getmechanism(:tippetop, Δt = 0.01, g = -9.00, contact = true);
initialize!(mech, :tippetop, x = [0,0,1.0], q = UnitQuaternion(RotX(0.99*π)), ω = [0,0.05,50.])
@elapsed storage = simulate!(mech, 5.0, record = true, verbose = false, opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)





function getrectangle(; Δt::T=0.01, g::T=-9.81) where T
    origin = Origin{T}(name="origin")
    mass = 1.0
    r = 0.1
    h = 1.0
    bodies = [Box(h/25, h/2, h, mass, color = RGBA(1., 0., 0.), name = "box")]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = "floating_joint")]
    mechanism = Mechanism(origin, bodies, eqcs, Δt = Δt, g = g)
    return mechanism
end

function initializerectangle!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}

    eqc = geteqconstraint(mechanism, "floating_joint")
    zeroVelocity!(mechanism)
    setPosition!(mechanism, eqc, [x; rotation_vector(q)])
    setVelocity!(mechanism, eqc, [v; ω])
end


mech = getmechanism(:rectangle, Δt = 0.01, g = -0.00);
body1 = collect(mech.bodies)[1]
body1.J = Matrix(Diagonal([1,2,3.]))
initialize!(mech, :rectangle, x = [0,0,1.0], q = UnitQuaternion(RotX(0.00*π)), ω = [0,5.0,0.01])
@elapsed storage = simulate!(mech, 8.0, record = true, verbose = false, opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)
