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


function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.joints)[2:end])
        nu = control_dimension(eqc)
        u = 50*mechanism.Δt*(ones(nu) .- 0.5)
        set_input!(eqc, u)
    end
    return
end

mech = getmechanism(:hopper, Δt = 0.05, g = -0*9.81, contact = false, limits = true,
    contact_body = true, spring = 1.0, damper = 0.0);
mech.ϕreg .= 0.0
initialize!(mech, :hopper, x = 0.0, z = 0.0, θ = -0.0)
@elapsed storage = simulate!(mech, 0.2, controller!, record = true,
    opts=InteriorPointOptions(verbose=true, btol = 1e-6))
visualize(mech, storage, vis = vis, show_contact = true)


set_entries!(mech)





plot(hcat(Vector.(storage.ω[1])...)')



eqc1 = collect(mech.joints)[1]
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]
zerodimstaticadjoint(∂g∂ʳpos(mech, eqc1, body1))
collect(mech.joints)



const Dojo = Main

# angular_damping!(mech, body1)
# ∂angular_damping!(mech, body1)

∇ = [szeros(3,6); sones(3,6)]
∇ = hcat(szeros(3,6), sones(3,6))
integrator_jacobian_velocity(one(UnitQuaternion), srand(3), 0.1)
#
#     for i=1:Nc
#         bnd = ineqc.constraints[i]
#         bnd_type = typeof(ineqc.constraints[i])
#
#         M = integrator_jacobian_velocity(q2, ω2, Δt)
#         function d(vars)
#             x = vars[1:3]
#             q = UnitQuaternion(vars[4:7]..., false)
#             return ∂g∂ʳpos(bnd, x, q, nothing)' * ineqc.γsol[2]
#         end
#
#         if bnd_type <: NonlinearContact
#             body.state.D -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
#         elseif bnd_type <: ImpactContact
#             body.state.D -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
#         elseif bnd_type <: LinearContact
#             body.state.D -= FiniteDiff.finite_difference_jacobian(d, [x3; q3.w; q3.x; q3.y; q3.z]) * M
#         end
#     end
#     return
# end


angular_damping!(mech, body1)

body1.state.q2
body2.state.q2


mech = getmechanism(:slider, Δt = 0.05, g = -9.81, spring = 10.0, damper = 1.0);
initialize!(mech, :slider, z1 = 0.0)
@elapsed storage = simulate!(mech, 3.00, record = true, verbose = false,
    opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis, show_contact = true)







# env = make("halfcheetah", vis = vis)

env.aspace
seed(env, s = 11)
obs = reset(env)[2]
render(env)

1000*sample(env.aspace)
collect(env.mechanism.joints)[1]
for i = 1:25
    render(env)
    sleep(0.05)
    # action = 120*env.mechanism.Δt*ones(6)#1000*sample(env.aspace) # your agent here (this takes random actions)
    action = sample(env.aspace)#1000*sample(env.aspace) # your agent here (this takes random actions)
    obs, r, done, info = step(env, action)
    @show r

    if done
        observation = reset(env)
    end
end
close(env)

env.mechanism.joints
control_dimension(env.mechanism)
sample(env.aspace)
# sample(env.aspace)
#
m.body_inertia
@show m.body_mass

# initialize!(env.mechanism, :halfcheetah, z = 2.0)
# torso = get_body(env.mechanism, "torso")
# eqc1 = get_joint_constraint(env.mechanism, "floating_joint")
# torso.state.x2
# orig = env.mechanism.origin
# minimal_coordinates(eqc1.constraints[1], orig, torso)
# minimal_coordinates(eqc1.constraints[2], orig, torso)


getMinState(env.mechanism)


env.x .= getMinState(env.mechanism)
render(env)

################################################################################
# Sparsify
################################################################################

using LinearAlgebra

nx = 5
nr = 10
nu = 5
Δt = 0.1
Rx0 = rand(nr, nx)
Ru0 = rand(nr, nu)
Rz1 = rand(nr, nr)
A = (Rz1 \ Rx0)[1:nx,:]
B = (Rz1 \ Ru0)[1:nx,:]

function idynamics(x1, x0, u0)
    return A*x0 + B*u0 - x1
end

function edynamics(x0, u0)
    return A*x0 + B*u0
end

x0 = rand(nx)
u0 = rand(nu)

x1 = edynamics(x0, u0)

M = [zeros(nr, nx+nu) inv(Rz1);
     Rx0 Ru0          1*Diagonal(ones(nr));
     ]
#    x0 u0            r0                    z1
M = [zeros(nr, nx+nu) zeros(nr, nr)         Diagonal(ones(nr)) ; # z1
     Rx0 Ru0          1*Diagonal(ones(nr))  Rz1                ; # r1
     ]

M
z1r1 = M \ [x0; u0; zeros(nr); z1]
z1 = z1r1[1:nr]
r1 = z1r1[nr .+ (1:nr)]
x1 = z1[1:nx]
norm(x1 - edynamics(x0, u0))

M = zeros(10,10)
for k = 1:10
    M[k,k] += rand()
end
for k = 1:9
    M[k+1,k] += rand()
    M[k,k+1] += rand()
end

M

inv(M)
