################################################################################
# Setup
################################################################################
# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

################################################################################
# Build mechanism and Identify A and B
################################################################################

mech = getmechanism(:atlas, timestep = 0.01, g = -9.81, cf = 0.8, contact = true)
initialize!(mech, :atlas, tran = [0,0,1.9291], rot = [0.,0,0])
for (i,joint) in enumerate(mech.joints)
    jt = joint.constraints[1]
    jr = joint.constraints[2]
    joint.isdamper = true #false
    joint.isspring = false #false

    jt.spring = 1 * 1.0 * 1e-0 .* sones(3)[1]# 1e4
    jt.damper = 1 * 1.0 * 1e-0 .* sones(3)[1]# 1e4
    jr.spring = 1 * 1.0 * 1e-0 .* sones(3)[1]# 1e4
    jr.damper = 1 * 1.0 * 1e+3 .* sones(3)[1]# 1e4
end

@elapsed storage = simulate!(mech, 0.1, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

# show sign distance function
contacts = collect(mech.contacts)
for (i,ineqc) in enumerate(contacts)
    ineqc = contacts[1]
    cont = ineqc.constraints[1]
    body = get_body(mech, ineqc.parentid)
    x3, q3 = current_configuration(body.state)
    sdf = cont.ainv3 * (x3 + vrotate(cont.p,q3) - cont.offset)
    println("sdf $i:", sdf)
end

#################################################################################
# get control matrices

# Set data
Nb = length(mech.bodies)
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
attjac = attitude_jacobian(data, Nb)

# IFT
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-6
plot(Gray.(abs.(1e10 .* datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(fd_solmat)))

solmat = full_matrix(mech.system)
indx2 = vcat([12i .+ [1,2,3,7,8,9] for i = 1:31]...)
indx3 = 30 * 5 .+ (1:31*6)
indu = 31*12 + 6 .+ (1:30)
A = - (solmat \ datamat[:,indx2])[indx3,:]
B = - (solmat \ datamat[:,indu])[indx3,:]
# size data 31 * 12 + 6 +  30 * 1
# size residual 30 * 5 + 31 * 6 + 8 * 8
Q = 10 * Matrix(Diagonal(ones(31*6)))
R = 1 * Matrix(Diagonal(ones(30)))
P = dare(A, B, Q, R)
K = R \ B' * P
cond(K)

plot(Gray.(abs.(K ./ 1e14)))

# PD control law
nu = sum([control_dimension(eqc, floatingbase = false) for eqc in collect(mech.joints)])
angles = [minimal_coordinates(mech, joint)[1] for joint in collect(mech.joints)[2:end]]
δangles = zeros(nu)
ind = 23
# δangles[ind] += π/2
angles += δangles

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.joints)[2:end])
        if control_dimension(joint) == 1
            θ = minimal_coordinates(mechanism, joint)[1]
            dθ = minimal_velocities(mechanism, joint)[1]
            u = 3e+2 * (angles[i] - θ) #+ 5e-2 * (0 - dθ)
            u = clamp(u, -150.0, 150.0) * mechanism.timestep
            if joint.name ∈ ("r_leg_akx", "r_leg_aky", "l_leg_akx", "l_leg_aky", "back_bkx", "back_bky", "back_bkz")
                u = 1e+2 * (angles[i] - θ) #+ 5e-2 * (0 - dθ)
                u = clamp(u, -100.0, 100.0) * mechanism.timestep
            end
            u = 0.0
            set_input!(joint, SA[u])
        end
    end
    return
end

# forcedstorage = simulate!(tmech, 2.5, controller!, record = true, solver = :mehrotra!)
# @elapsed forcedstorage = simulate!(tmech, 2.5, controller!, record = true, solver = :mehrotra!)
# @elapsed forcedstorage = simulate!(mech, 1.5, controller!, record = true, solver = :mehrotra!)
# @profiler forcedstorage = simulate!(tmech, 0.5, controller!, record = true, solver = :mehrotra!)
# visualize(tmech, forcedstorage, vis = vis)

@elapsed storage = simulate!(mech, 4, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)



gains = zeros(30, 2)
gains[23,:] = [1e-1, 5e-2]

nams = [eqc.name for eqc in mech.joints]

nams[1:10]
nams[11:20]
nams[21:30]

# Set data
Nb = length(mech.bodies)
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
attjac = attitude_jacobian(data, Nb)

# IFT
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-6
plot(Gray.(abs.(1e10 .* datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵ = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 8e-3
plot(Gray.(1e10 .* sensi))
plot(Gray.(fd_sensi))



function dare(A, B, Q, R)
    n = size(A, 1);

    E = [
        Matrix{Float64}(I, n, n) B*(R\B');
        zeros(size(A)) A'
    ];
    F = [
        A zeros(size(A));
        -Q Matrix{Float64}(I, n, n)
    ];

    QZ = schur(F, E);
    QZ = ordschur(QZ, abs.(QZ.alpha./QZ.beta) .< 1);

    return QZ.Z[(n+1):end, 1:n]/QZ.Z[1:n, 1:n];
end
