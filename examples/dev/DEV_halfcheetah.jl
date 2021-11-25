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


function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = controldim(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.eqconstraints
        nu = controldim(eqc)
        if eqc.parentid != nothing
            body = getbody(mechanism, eqc.parentid)
            rot = eqc.constraints[2]
            A = Matrix(nullspacemat(rot))
            Fτ = springforce(mechanism, eqc, body)
            F = Fτ[1:3]
            τ = Fτ[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end

function controller!(mechanism, k,)
    u = [zeros(3); 0.4; zeros(5)]
    # u = [zeros(4); 0.4; zeros(5)]
    # u = [zeros(5); 0.4; zeros(5)]
    off = 0
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
        nu = controldim(eqc)
        setForce!(mechanism, eqc, SVector{nu}(u[off .+ (1:nu)]))
        off += nu
    end
    return
end

mech = getmechanism(:halfcheetah, Δt = 0.01, g = -9.81, damper = 1.0, spring = 00.0, contact = true)
initialize!(mech, :halfcheetah, x = 0.0, z = 0.0, θ = -0.0)


@elapsed storage = simulate!(mech, 10.05, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
get_sdf(mech, storage)[1][1]
get_sdf(mech, storage)[2][1]




ugc = gravity_compensation(mech)

function controller!(mechanism, k,)
    u = ugc
    off = 0
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
        nu = controldim(eqc)
        setForce!(mechanism, eqc, SVector{nu}(u[off .+ (1:nu)]))
        off += nu
    end
    return
end

getfield.(mech.bodies, :m)


mech = getmechanism(:halfcheetah, Δt = 0.01, g = -9.81, damper = 10.0, spring = 0.0)
initialize!(mech, :halfcheetah, x = 0.0, z = 0.00, θ = 0.0)
@elapsed storage = simulate!(mech, 7.0, controller!, record = true, solver = :mehrotra!, verbose = false)

visualize(mech, storage, vis = vis)
