"""
    Total mechanical energy of a mechanism.
"""
function mechanicalEnergy(mechanism::Mechanism)
    ET = kineticEnergy(mechanism)
    EV = potentialEnergy(mechanism)
    E = ET + EV
    return E
end

"""
    Kinetic energy of a mechanism due to translation and rotation velocities.
"""
function kineticEnergy(mechanism::Mechanism)
    Δt = mechanism.Δt
    ET = 0.0
    for body in mechanism.bodies
        x2, v2, q2, ω2 = fullargssol(body.state)
        p1 = [body.m * v2; body.J * ω2] # linear momentum, angular momentum
        p1 += 0.5 * [0, 0, body.m * mechanism.g * Δt, 0, 0, 0] # gravity
        for eqc in mechanism.eqconstraints
            if body.id == eqc.parentid
                for (i,joint) in enumerate(eqc.constraints)
                    cbody = getbody(mechanism, eqc.childids[i])
                    p1 += 0.5 * springforcea(joint, body.state, cbody.state, Δt)
                end
            end
            for (i,joint) in enumerate(eqc.constraints)
                if eqc.childids[i] == body.id
                    if eqc.parentid != nothing
                        pbody = getbody(mechanism, eqc.parentid)
                        p1 += 0.5 * springforceb(joint, pbody.state, body.state, Δt)
                    else
                        p1 += 0.5 * springforceb(joint, body.state, Δt)
                    end
                end
            end
        end
        M = [body.m * I szeros(3,3); szeros(3,3) body.J]
        v1 = M \ p1
        ET += 0.5 * v1'* M * v1
    end
    return ET
end


function oldkineticEnergy(mechanism::Mechanism)
    ET = 0.0
    for body in mechanism.bodies
        spring = collect(mechanism.eqconstraints)[1].constraints[1].spring
        ET += kineticEnergy(body, spring, mechanism.g, mechanism.Δt)
    end
    return ET
end

"""
    Kinetic energy of a body due to translation and rotation velocities.
"""
function oldkineticEnergy(body::Body, spring, g, Δt)
    # x1, v1, q1, ω1 = fullargsc(body.state) # TODO we should compute the velocities at t = 2 not t = 2.5
    # ET1 = 0.5 * body.m * v1'*v1
    # ET1 += 0.5 * ω1' * body.J * ω1
    # x2, v2, q2, ω2 = fullargssol(body.state) # TODO we should compute the velocities at t = 2 not t = 2.5
    # ET2 = 0.5 * body.m * v2'*v2
    # ET2 += 0.5 * ω2' * body.J * ω2
    x2, v2, q2, ω2 = fullargssol(body.state) # TODO we should compute the velocities at t = 2 not t = 2.5
    # p1 = body.m * v2 + Δt * [0, 0, body.m * g]
    x2, q2 = posargsk(body.state)
    z = x2[3]
    p1 = body.m * v2
    p1 += - Δt/2 * spring * [0, 0, z + 0.5]
    p1 += - Δt/2 * [0, 0, -body.m * g]
    v1 = p1 / body.m
    ET = 0.5 * body.m * v1'*v1
    # # ET += 0.5 * ω1' * body.J * ω1
    # return (ET1 + ET2) / 2
    # return ET2
    return ET
end


"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potentialEnergy(mechanism::Mechanism)
    V = 0.0
    for body in mechanism.bodies
        V += potentialEnergy(mechanism, body)
    end
    for eqc in mechanism.eqconstraints
        V += potentialEnergy(mechanism, eqc)
    end
    return V
end

"""
    Potential energy of one body due to gravity.
"""
function potentialEnergy(mechanism::Mechanism, body::Body)
    x, q = posargsk(body.state)
    z = x[3]
    V = - body.m * mechanism.g * z
    return V
end

"""
    Potential energy of a joint due to the spring.
"""
function potentialEnergy(mechanism::Mechanism, eqc::EqualityConstraint)
    V = 0.0
    parentid = eqc.parentid
    Δt = mechanism.Δt
    if parentid != nothing
        pbody = getbody(mechanism, parentid)
    else
        pbody = mechanism.origin
    end
    for (i,joint) in enumerate(eqc.constraints)
        cbody = getbody(mechanism, eqc.childids[i])
        V += potentialEnergy(joint, pbody, cbody, Δt)
    end
    return V
end

"""
    Potential energy of a joint due to the spring.
"""
potentialEnergy(joint::Rotational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Rotational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Origin, cbody::Body, Δt::T) where {T} = 0.0
potentialEnergy(joint::Translational{T,3}, pbody::Body, cbody::Body, Δt::T) where {T} = 0.0

function potentialEnergy(joint::Translational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Translational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(1,2,3)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Body, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, pbody.state, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end

function potentialEnergy(joint::Rotational, pbody::Origin, cbody::Body, Δt::T) where {T}
    if joint.spring > 0.0
        force = springforceb(joint, cbody.state, Δt)[SVector{3,Int}(4,5,6)] / Δt
        return 0.5 * force'*force ./ joint.spring
    elseif joint.spring == 0.0
        return 0.0
    else
        @warn "invalid potential energy: negative spring constant."
        return 0.0
    end
end


vis = Visualizer()
open(vis)

# mech = getmechanism(:snake, contact = false, spring = 000.0)
# initialize!(mech, :snake)
# storage = simulate!(mech, 10.1, record = true, solver = :mehrotra!, verbose = false)
# visualize(mech, storage, vis = vis)


# body1 = collect(mech.bodies)[1]
# eqc1 = collect(mech.eqconstraints)[1]
# eqc2 = collect(mech.eqconstraints)[2]
# eqc1.parentid
# eqc2.parentid
# eqc2.constraints[1].spring

potentialEnergy(mech)
kineticEnergy(mech)
mechanicalEnergy(mech)

function get_energy(tsim, Δt)
    # mech = getmechanism(:slider, Δt = Δt, g = -1.0, spring = 100.0)
    # initialize!(mech, :slider, z1 = 2.0)
    mech = getmechanism(:nslider, Δt = Δt, g = -9.81, spring = 10.0, Nlink = 5)
    initialize!(mech, :nslider, Δz = 1.0)
    telap = @elapsed simulate!(mech, tsim, record = false, solver = :mehrotra!, verbose = false)

    ET = kineticEnergy(mech)
    EV = potentialEnergy(mech)
    EM = mechanicalEnergy(mech)
    return [ET, EV, EM, tsim / telap]
end

ts = [0.11 * i for i = 1:50]
energies = [get_energy(t, 0.02) for t in ts]
plt = plot()
plot!(ts, hcat(energies...)[1,:], linewidth = 3.0, label = "kinetic")
plot!(ts, hcat(energies...)[2,:], linewidth = 3.0, label = "potential")
plot!(ts, hcat(energies...)[3,:], linewidth = 3.0, label = "mechanical")
plt = plot()
energies[1][3]
err = abs.((hcat(energies...)[3,:] .- energies[1][3]) / energies[1][3])
plot(ts, log.(10, err), linewidth = 3.0, label = "mechanical", ylims = (-8, 2))


Δts = [1e-6*exp(1/4 * log(10))^(i-1) for i = 1:20]
energies = [get_energy(Δt*897, Δt) for Δt in Δts]
plt = plot()
err = abs.((hcat(energies...)[3,:] .- energies[10][3]) / energies[10][3])
tratio = hcat(energies...)[4,:]
scatter(log.(10, tratio), log.(10, err), linewidth = 3.0, label = "mechanical", ylims = (-10, 2))

mech = getmechanism(:pendulum)
initialize!(mech, :pendulum, ϕ1 = π/2)
storage = simulate!(mech, 2.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


mech = getmechanism(:slider, Δt = 0.02, g = 0.0, spring = 2.0)
initialize!(mech, :slider, z1 = 1)

mech = getmechanism(:nslider, Δt = 0.01, g = -10.0, spring = 100.0, Nlink = 5)
initialize!(mech, :nslider, z1 = 0.0, Δz = 0.5)
storage = simulate!(mech, 2.5, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


body1 = collect(mech.bodies)[1]
body1.state
eqc1 = collect(mech.eqconstraints)[1]
tra1 = eqc1.constraints[1]
tra1.spring
potentialEnergy(mech, eqc1)
0.5 * 100 * z*2

gv = mech.g
m = body1.m
x, q = posargsk(body1.state)
z = x[3]
EV = potentialEnergy(mech, body1)
potentialEnergy(mech, body1) + m * gv * z
