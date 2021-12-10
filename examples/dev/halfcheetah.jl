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


mech = getmechanism(:halfcheetah, Δt = 0.10, g = -9.81, contact = true);
initialize!(mech, :halfcheetah, x = 0.0, z = 0.0, θ = -0.0)
@elapsed storage = simulate!(mech, 0.30, record = true, verbose = false,
    opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)

env = make("halfcheetah")
open(env.vis)

seed(env, s = 11)
obs = reset(env)[2]
render(env)


collect(env.mechanism.eqconstraints)[1]
for i = 1:10
    render(env)
    action = 1000*sample(env.aspace) # your agent here (this takes random actions)
    obs, r, done, info = step(env, action)
    @show r

    if done
        observation = reset(env)
    end
end
close(env)


# sample(env.aspace)
#


# initialize!(env.mechanism, :halfcheetah, z = 2.0)
# torso = getbody(env.mechanism, "torso")
# eqc1 = geteqconstraint(env.mechanism, "floating_joint")
# torso.state.x2
# orig = env.mechanism.origin
# minimalCoordinates(eqc1.constraints[1], orig, torso)
# minimalCoordinates(eqc1.constraints[2], orig, torso)


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
