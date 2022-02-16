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


mech = getmechanism(:ant, timestep = 0.01, g = -9.81, contact = true,
    contact_body = true, spring = 0.0, damper = 1.0);
initialize!(mech, :ant, rot = [0,0,0.], ankle = 0.25)
@elapsed storage = simulate!(mech, 2.0, record = true, verbose = false,
    opts=SolverOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)

env = make("ant", vis = vis)

env.aspace
seed(env, s = 11)
obs = reset(env)[2]
render(env)

1000*sample(env.aspace)
collect(env.mechanism.joints)[1]
for i = 1:25
    render(env)
    sleep(0.05)
    # action = 120*env.mechanism.timestep*ones(6)#1000*sample(env.aspace) # your agent here (this takes random actions)
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








using BenchmarkTools
eqc1 = collect(mech.joints)[1]
eqc2 = collect(mech.joints)[2]
joint21 = eqc2.constraints[1]
joint22 = eqc2.constraints[2]
@benchmark joint_impulse_index(eqc1, 1)
@benchmark joint_impulse_index($eqc2, $1)


@inline function correction(joint::Joint{T,Nλ,Nb,N}, Δ, μ) where {T,Nλ,Nb,N}
    Δs, Δγ = split_impulses(joint, Δ)
	return [- Δs .* Δγ .+ μ; szeros(Nb + Nλ)]
end


# μ = mech.μ
# jj = fill(joint22, 2)
# ll = fill(λ, 2)
# @benchmark vv($μ, $ll, $jj)
# @code_warntype vv(μ, ll, jj)

λ = eqc2.λsol[2]
μ = mech.μ

system = mech.system
residual_entries = mech.residual_entries
res = residual_entries[2]
ste = get_entry(system, 2)
@benchmark correction!($mech, $res, $ste, $eqc2)




correction!(mech, λ, λ, eqc2)

jj = eqc2.constraints[]1
ll = fill(λ,2)
@benchmark correction.($jj, $ll, $μ)
@benchmark correction($joint22, $λ, $μ)

function ttt()
	c = szeros(0)
	for i = 1:10
		c = vcat(c, sones(10))
	end
	return c
end

@benchmark ttt()
@code_warntype ttt()




system = mech.system
residual_entries = mech.residual_entries
res = residual_entries[2]
ste = get_entry(system, 2)
@benchmark correction($mech, $ste, $eqc2)
@code_warntype correction(mech, ste, eqc2)

@benchmark correction!($mech, $res, $ste, $eqc2)
@code_warntype correction!(mech, res, ste, eqc2)


@benchmark ∂g∂ʳself($joint22, $λ)
@code_warntype ∂g∂ʳself(joint22, λ)

joint22
SMatrix{2,2,Float64,4}([1 2; 3 4])




@benchmark ∂g∂ʳself($joint22, $λ)
∂g∂ʳself(joint22, λ)
