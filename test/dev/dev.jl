
vis = Visualizer()
open(vis)

mech = getmechanism(:snake, contact = false, spring = 000.0, g = 0.0)
initialize!(mech, :snake, v = zeros(3), Δv = zeros(3), Δω = zeros(3))
storage = simulate!(mech, 10.1, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

include("conservation_test.jl")

potentialEnergy(mech)
kineticEnergy(mech)
mechanicalEnergy(mech)

function get_energy(tsim, Δt)
    # mech = getmechanism(:slider, Δt = Δt, g = -1.0, spring = 100.0)
    # initialize!(mech, :slider, z1 = 2.0)
    mech = getmechanism(:nslider, Δt = Δt, g = -9.81, spring = 10.0, Nb = 5)
    initialize!(mech, :nslider, Δz = 1.0)
    telap = @elapsed simulate!(mech, tsim, record = false, solver = :mehrotra!, verbose = false)

    ET = kineticEnergy(mech)
    EV = potentialEnergy(mech)
    EM = mechanicalEnergy(mech)
    return [ET, EV, EM, tsim / telap]
end

ts = [0.10 * i for i = 1:50]
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


mech = getmechanism(:nslider, Δt = 0.01, g = -10.0, spring = 100.0, Nb = 5)
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
x, q = posargs2(body1.state)
z = x[3]
EV = potentialEnergy(mech, body1)
potentialEnergy(mech, body1) + m * gv * z


################################################################################
# atlas
################################################################################

include("conservation_test.jl")
Random.seed!(100)
Δt_ = 0.01
mech = getmechanism(:atlas, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.05, contact = false)
initialize!(mech, :atlas, tran = [0, 0, 2.])

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if controldim(joint) == 1
            if k ∈ (1:50)
                u = 1e0 * (1.0 - 0.5) * Δt_
            elseif k ∈ (51:100)
                u = -1e0 * (1.0 - 0.5) * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

# mech.eqconstraints[1].constraints[1].spring
# mech.eqconstraints[1].constraints[2].spring
# mech.eqconstraints[1].constraints[1].damper
# mech.eqconstraints[1].constraints[2].damper
momentum(mech)

[eqc.name for eqc in mech.eqconstraints]
[(body.name, body.m) for body in mech.bodies]

################################################################################
# snake
################################################################################
include("conservation_test.jl")
n = 1000
Δt_ = 0.1 / n
Nb_ = 2
mech = getmechanism(:snake, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0, contact = false, Nb = Nb_)
initialize!(mech, :snake, Δv = zeros(3), Δω = zeros(3))

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if  i == Nb_-1
            if controldim(joint) == 1
                # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
                if k ∈ (1:n)
                    u = 1e-0 * Δt_
                elseif k ∈ (n+1:2n)
                    u = -1e-0 * Δt_
                else
                    u = 0.0
                end
                setForce!(mechanism, joint, SA[u])
            end
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)

momentum(mech, collect(mech.bodies)[1])[4]
momentum(mech, collect(mech.bodies)[2])[4]
momentum(mech)[4]

plot([ω[1] for ω in storage.ω[1]])
plot([ω[1] for ω in storage.ω[2]])

collect(mech.bodies)[1].m

bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)
body1 = bodies[1]
eqc1 = eqcs[1]
tra1 = eqc1.constraints[1]
rot1 = eqc1.constraints[2]
eqc2 = eqcs[2]
tra2 = eqc2.constraints[1]
rot2 = eqc2.constraints[2]

rot1.qoffset
rot2.qoffset


################################################################################
# slider
################################################################################
include("conservation_test.jl")
Δt_ = 0.001
mech = getmechanism(:slider, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0)
initialize!(mech, :slider)

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[1:end])
        if controldim(joint) == 1
            # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
            if k ∈ (1:500)
                u = 5e-1 * Δt_
            elseif k ∈ (501:1000)
                u = -5e-1 * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

momentum(mech)

plot([x[3] for x in storage.x[1]])
plot([ω[1] for ω in storage.ω[1]])


################################################################################
# npendulum
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nb_ = 2
mech = getmechanism(:npendulum, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0, Nb = Nb_)
initialize!(mech, :npendulum, ϕ1 = 0.)#, Δv = zeros(3), Δω = zeros(3))

function controller!(mechanism, k)
    n = 50
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if  i == Nb_-1
            if controldim(joint) == 1
                # u = 5e-1 * (rand() - 0.5) * Δt_ * (k < 100)
                if k ∈ (1:n)
                    u = 1e1 * Δt_
                elseif k ∈ (n+1:2n)
                    u = -1e1 * Δt_
                else
                    u = 0.0
                end
                setForce!(mechanism, joint, SA[u])
            end
        end
    end
    return
end

qq = UnitQuaternion(rand(4)...)
τ0 = rand(3)
τ1 = vrotate(vrotate(τ0, qq), inv(qq))
τ0 - τ1


# storage = simulate!(mech, 50*2Δt_, controller!, record = true, solver = :mehrotra!, verbose = true)
storage = simulate!(mech, 2.0, controller!, record = true, solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)

momentum(mech, collect(mech.bodies)[1])[4]
momentum(mech, collect(mech.bodies)[2])[4]

plot([ω[1] for ω in storage.ω[1]])
plot([ω[1] for ω in storage.ω[2]])

bodies = collect(mech.bodies)
eqcs = collect(mech.eqconstraints)
body1 = bodies[1]
eqc1 = eqcs[1]
tra1 = eqc1.constraints[1]
rot1 = eqc1.constraints[2]
eqc2 = eqcs[2]
tra2 = eqc2.constraints[1]
rot2 = eqc2.constraints[2]

rot1.qoffset
rot2.qoffset




vis = Visualizer()
open(vis)





################################################################################
# humanoid
################################################################################

include("conservation_test.jl")
Random.seed!(100)
Δt_ = 0.01
mech = getmechanism(:humanoid, Δt = Δt_, g = 0.0, spring = 0.0, damper = 3., contact = false)
initialize!(mech, :humanoid, tran = [0, 0, 2.])

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if controldim(joint) == 1
            if k ∈ (50:150)
                u = 2e0 * (1.0 - 0.5) * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 100.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, downsample(storage, 100), vis = vis)

# storage2500 = deepcopy(storage)
# storage500 = deepcopy(storage)

norm(momentum(mech)[4:6])
momentum(mech)



################################################################################
# snake
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nb_ = 2

Random.seed!(100)
ω_ = 0.0*rand(3)
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 0.0*rand(3)

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = controldim(joint)
        if nu <= 5
            if k ∈ (50:150)
                u = 2e0 * (1.0 - 0.5) * Δt_ * [1.0, 0.0, 1.0] #[0.0; 1.0; zeros(nu-2)]
            else
                u = zeros(nu)
            end
            setForce!(mechanism, joint, SA[u...])
        end
    end
    return
end
mech.eqconstraints[1]
mech.eqconstraints[2]
mech = getmechanism(:snake, Δt = Δt_, g = 0.0, spring = 0.0, damper = 1.0, contact = false, Nb = Nb_, jointtype = :Spherical)
initialize!(mech, :snake, ω = ω_, v = v_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 10.0, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)


visualize(mech, downsample(storage, 1), vis = vis)


################################################################################
# snake initial velocity
################################################################################
vis = Visualizer()
open(vis)

include("conservation_test.jl")
n = 1
Δt_ = 0.01
Δt_ /= n
Nb_ = 2

Random.seed!(100)
ω_ = -1e-1*[1,1,1.0] #* 100 * Δt_
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 0.0*4.0*[0,0,1.0] #* 100 * Δt_

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = controldim(joint)
        if nu <= 5
            if k ∈ (10:10 + 1000n)
                u = 1.0 * 3e-1 * Δt_ * [0.1, 1.0, 1.0] #[0.0; 1.0; zeros(nu-2)]
            elseif k ∈ (10 + 1000n:10 + 2000n)
                u = 0.0 * -3e-1 * Δt_ * [1.0, 0.0, 1.0] #[0.0; 1.0; zeros(nu-2)]
            else
                u = zeros(3)
            end
            setForce!(mechanism, joint, SA[u...])
        end
    end
    return
end
mech = getmechanism(:snake, Δt = Δt_, g = -0.05, spring = 0.0, damper = 1.00, contact = false, Nb = Nb_, jointtype = :Spherical)
initialize!(mech, :snake, ω = ω_, v = v_, Δv = Δv_, Δω = Δω_, ϕ1 = 0.0)
storage = simulate!(mech, 100.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)
visualize(mech, downsample(storage, 10), vis = vis)

function get_momentum(h)
    mech = getmechanism(:snake, Δt = Δt_, g = -0.05, spring = 0.0, damper = 1., contact = false, Nb = Nb_, jointtype = :Spherical)
    initialize!(mech, :snake, ω = ω_, v = v_, Δv = Δv_, Δω = Δω_, ϕ1 = 0.0)
    storage = simulate!(mech, h, controller!, record = true, solver = :mehrotra!, verbose = false)
    m0 = momentum(mech)
    return m0[4:6]
end
mm = [get_momentum(0.1+0.5i) for i = 1:2]
tt = [0.1+0.5i for i = 1:2]
mm = [i .- mm[1] for i in mm]
plot(tt, hcat(mm...)')


eqc2 = collect(mech.eqconstraints)[2]
tra2 = eqc2.constraints[1]
tra2.vertices

eqc = mech.eqconstraints[2]
f1 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[3])) * eqc.λsol[2])[1:3]
f2 = (zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[4])) * eqc.λsol[2])[1:3]
norm(f1 + f2)
t1 = vrotate((zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[3])) * eqc.λsol[2])[4:6], mech.bodies[3].state.q2[1])
t2 = vrotate((zerodimstaticadjoint(∂g∂ʳpos(mech, eqc, mech.bodies[4])) * eqc.λsol[2])[4:6], mech.bodies[4].state.q2[1])
norm(t1 + t2)
t1
t2

bodya.J
bodya = collect(mech.bodies)[1]
bodyb = collect(mech.bodies)[2]
xa, qa = posargs2(bodya.state)
xb, qb = posargs2(bodyb.state)
rotation_matrix(qa) * tra2.vertices[1]
rotation_matrix(qb) * tra2.vertices[2]
tra2.vertices[1]
tra2.vertices[2]
qa
qb

rotation_matrix(qa) * (10*I) * rotation_matrix(inv(qa))

################################################################################
# pendulum
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nb_ = 2

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if controldim(joint) == 1
            if k ∈ (1:10)
                u = 1e0 * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

Random.seed!(100)
ω_ = 0.0*rand(3)
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 0.0*rand(3)
ϕ1_ = 0.0
mech = getmechanism(:npendulum, Δt = Δt_, g = 0.00, spring = 0.0, damper = 0.0, Nb = Nb_)
initialize!(mech, :npendulum, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 0.50, controller!, record = true, solver = :mehrotra!, verbose = true)
m0 = momentum(mech)

mech = getmechanism(:npendulum, Δt = Δt_, g = 0.00, spring = 0.0, damper = 0.0, Nb = Nb_)
initialize!(mech, :npendulum, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 1.00, controller!, record = true, solver = :mehrotra!, verbose = true)
m1 = momentum(mech)
norm((m1 - m0)[4:6], Inf)
(m1 - m0)[5]
m1 - m0
visualize(mech, storage, vis = vis)

eqc2 = collect(mech.eqconstraints)[2]
body2 = collect(mech.bodies)[2]
zerodimstaticadjoint(∂g∂ʳpos(mech, eqc2, body2)) * eqc2.λsol[2]
∂g∂ʳpos(mech, eqc2, body2)
eqc2.λsol[2]



################################################################################
# snake
################################################################################
include("conservation_test.jl")
Δt_ = 0.01
Nb_ = 2

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = controldim(joint)
        if 5 >= nu >= 1
            if k ∈ (1:40)
                u = 3e-0 * Δt_ * [1.0; zeros(nu-1)]
            else
                u = 0.0 * [1.0; zeros(nu-1)]
            end
            setForce!(mechanism, joint, SA[u...])
        end
    end
    return
end

Random.seed!(100)
ω_ = 0.0*rand(3)
v_ = 0.0*rand(3)
Δv_ = 0.0*rand(3)
Δω_ = 0.0*rand(3)
ϕ1_ = pi/2
jointtype = :Fixed
mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.3, Nb = Nb_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 1.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m0 = momentum(mech)

mech = getmechanism(:snake, Δt = Δt_, g = 0.00, contact = false, spring = 0.0, damper = 0.3, Nb = Nb_, jointtype = jointtype)
initialize!(mech, :snake, ϕ1 = ϕ1_, Δv = Δv_, Δω = Δω_)
storage = simulate!(mech, 60.00, controller!, record = true, solver = :mehrotra!, verbose = false)
m1 = momentum(mech)
norm((m1 - m0)[4:6], Inf)
(m1 - m0)[5]
m1 - m0
visualize(mech, storage, vis = vis)




vis = Visualizer()
open(vis)




function getwalker2d(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = "examples/examples_files/walker2d.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=true, g = g)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints[2:end]))
        eqc.isdamper = true
        eqc.isspring = true
        for joint in eqc.constraints
            joint.spring = spring
            joint.damper = damper
        end
    end

    # if contact
    #     origin = Origin{T}()
    #     bodies = Vector{Body{T}}(collect(mech.bodies))
    #     eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))
    #
    #     # Foot contact
    #     contacts = [
    #         [-0.1; -0.05; 0.0],
    #         [+0.1; -0.05; 0.0],
    #         [-0.1; +0.05; 0.0],
    #         [+0.1; +0.05; 0.0],
    #         ]
    #     n = length(contacts)
    #     normal = [[0;0;1.0] for i = 1:n]
    #     cf = cf * ones(T, n)
    #
    #     contineqcs1 = contactconstraint(getbody(mech, "left_foot"), normal, cf, p = contacts)
    #     contineqcs2 = contactconstraint(getbody(mech, "right_foot"), normal, cf, p = contacts)
    #
    #     setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
    #     mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    # end
    return mech
end

# function initializewalker2d!(mechanism::Mechanism; x::T = 0.0, z::T = 0.0, θ::T = 0.0) where {T}
#     setPosition!(mechanism,
#                  geteqconstraint(mechanism, "base_joint"),
#                  [x, z, θ])
# end
function initializewalker2d!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,1.2],
        rot::AbstractVector{T} = [0.1,0,0]) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "auto_generated_floating_joint"),
                 [tran; rot])
end
################################################################################
# walker2d
################################################################################

include("conservation_test.jl")
Random.seed!(100)
Δt_ = 0.01
mech = getmechanism(:walker2d, Δt = Δt_, g = 0.0, spring = 0.0, damper = 0.0, contact = false)
# initialize!(mech, :walker2d, x = 0.0, z = 0.0, θ = 0.0)
initialize!(mech, :walker2d, tran = [0, 0, 1.25], rot = zeros(3))
getfield.(collect(mech.eqconstraints), :name)
function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        if controldim(joint) == 1
            if k ∈ (1:100)
                u = 1e0 * (1.0 - 0.5) * Δt_
            else
                u = 0.0
            end
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

storage = simulate!(mech, 1.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
getfield.(collect(mech.bodies), :name)
