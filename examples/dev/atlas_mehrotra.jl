# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)


# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

timestep0 = 0.1
mech = getmechanism(:atlas, timestep = timestep0, g = -9.81, friction_coefficient = 0.8, contact = true,
    spring = 0.0, damper = 30.0, model_type = :fast)
initialize!(mech, :atlas, tran = [0,0,1.1], rot = [0.1,0.05,0])

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.joints)[2:end])
        pbody = get_body(mech, eqc.parentid)
        minJ = minimum(diag(pbody.J))
        for (i,joint) in enumerate(eqc.constraints)
            cbody = get_body(mech, eqc.childids[i])
            minJ = min(minJ, minimum(diag(cbody.J)))
        end
        nu = input_dimension(eqc)
        u = 1 * minJ * (rand(nu) .- 0.2) * timestep0 * 0.0
        set_input!(eqc, SVector{nu}(u))
    end
    return
end

@elapsed storage = simulate!(mech, 2.71, controller!, record = true, undercut = Inf,
    solver = :mehrotra!, verbose = true)
visualize(mech, storage, vis = vis)



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
plot(Gray.(abs.(1e10 .* fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(1e10 * fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 8e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))


test_solmat(:atlas, ϵ = 1e-8)
test_datamat(:atlas, ϵ = 1e-6)
test_sensitivity(:atlas, ϵ = 8e-3)

test_solmat(:quadruped, ϵ = 1e-8)
test_datamat(:quadruped, ϵ = 1e-6)
test_sensitivity(:quadruped, ϵ = 8e-3)
