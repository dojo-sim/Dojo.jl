function getpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
    spring = 0.0, damper = 0.0, spring_offset = szeros(1)) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1,
        joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset,
        rot_joint_limits = [SVector{1}([-0.25 * π]), SVector{1}([0.25 * π])]))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end

mech = getpendulum(Δt = 0.01, g = -9.81, spring = 0.0, damper = 0.0)
length(mech.eqconstraints[1].constraints[1])

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

# mech = getmechanism(:pendulum, Δt = 0.01, g = -9.81, spring = 100.0, damper = 5.0)
Random.seed!(100)
ϕ1 = 0.0 * π
initialize!(mech, :pendulum, ϕ1 = ϕ1)

function cont!(mechanism, k; u = 30.1)
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        nu = controldim(eqc, ignore_floating_base = false)
        su = mechanism.Δt * u * sones(nu)
        setForce!(mechanism, eqc, su)
    end
    return
end

∂g∂ʳvel(mech, [eq for eq in mech.eqconstraints][1], mech.bodies[2])
∂gab∂ʳba(mech, mech.bodies[2], [eq for eq in mech.eqconstraints][1])
∂g∂ʳposb(mech, [eq for eq in mech.eqconstraints][1], mech.bodies[2])
∂gab∂ʳba(mech, [eq for eq in mech.eqconstraints][1], mech.bodies[2])

function _∂g∂ʳposa(mechanism, eqc::EqualityConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    i = 2
    @show "hi"
    @show "yo yo"
    @show eqc.constraints[i]
    ∂g∂ʳposa(eqc.constraints[i], body, getbody(mechanism, eqc.childids[i]), eqc.childids[i], eqc.λsol[2][λindex(eqc,i)], mechanism.Δt)
end


_∂g∂ʳposa(mech, [eq for eq in mech.eqconstraints][1], mech.bodies[2])

joint_limits_length(mech.eqconstraints[1].constraints[2])



storage = simulate!(mech, 1.0, record = true, verbose = false)
visualize(mech, storage, vis = vis)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)

data = getdata(mech)
v15 = data[4:6]
sol = getsolution(mech)
v25 = sol[6:8]
norm(v15 - v25)

# IFT
setentries!(mech)
datamat = full_data_matrix(mech, attjac = true)
datamat1 = full_data_matrix(mech, attjac = false)
datamat2 = full_data_matrix(mech, attjac = false) * attjac
plot(Gray.(attjac))

solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show cond(solmat)
@show rank(solmat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
fd_datamat1 = finitediff_data_matrix(mech, data, sol)
@test norm(fd_datamat + datamat, Inf) < 1e-7
@test norm(fd_datamat1 + datamat1, Inf) < 1e-7
@test norm(fd_datamat + datamat2, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))


fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(fd_solmat)))
plot(Gray.(abs.(solmat)))


fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
norm(fd_sensi - sensi, Inf)
norm(fd_sensi, Inf)

###############################################################################
# plot
###############################################################################

plot(hcat(Vector.(storage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')


sdf = get_sdf(mech, storage)
plot(hcat(sdf[1]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[2]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[3]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[4]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[5]...)')
plot(hcat(sdf[6]...)')
plot(hcat(sdf[7]...)')
plot(hcat(sdf[8]...)')
