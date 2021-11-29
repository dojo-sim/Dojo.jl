jointtypes = [
    :Fixed,
    :Prismatic,
    :Planar,
    :FixedOrientation,
    :Revolute,
    :Cylindrical,
    :PlanarAxis,
    :FreeRevolute,
    :Orbital,
    :PrismaticOrbital,
    :PlanarOrbital,
    :FreeOrbital,
    :Spherical,
    :CylindricalFree,
    :PlanarFree
    ]
vis = Visualizer()
open(vis)


function vizMax(mech, z)
	Nb = length(mech.bodies)
	storage = Storage(2,Nb)
	for t = 1:2
		for b = 1:Nb
			x2, v15, q2, ϕ15 = unpackMaxState(z, b)
			storage.x[b][t] = x2
			storage.v[b][t] = v15
			storage.q[b][t] = q2
			storage.ω[b][t] = ϕ15
		end
	end
	visualize(mech, storage, vis = vis)
end


Δt0 = 0.01
g0 = -9.81
spring0 = 0.0
damper0 = 1.0
mech = getmechanism(:hopper, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :hopper, v = [0.1, 0.2, 0.3], ω = [0.4, 0.4, 0.1])

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
        nu = controldim(eqc)
        setForce!(mechanism, eqc, 0*SVector{nu}(Δt0 * (rand(nu) .- 0.5) ))
    end
    return
end
@elapsed storage = simulate!(mech, 3.00, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = 1e-14)
visualize(mech, storage, vis = vis)


################################################################################
# Test min2max max2min
################################################################################

function quaterror(q0::AbstractVector, q1::AbstractVector, p::Real = Inf)
	q0_ = UnitQuaternion(q0...)
	q1_ = UnitQuaternion(q1...)
	q_ = q0_ * inv(q1_)
	return norm([q_.x, q_.y, q_.z], p)
end

################################################################################
# Test setPosition! and setVelocity!
################################################################################
eqc1 = collect(mech.eqconstraints)[1]
eqc2 = collect(mech.eqconstraints)[2]
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]
tra2 = eqc2.constraints[1]
rot2 = eqc2.constraints[1]

x = srand(1)
Δx = zerodimstaticadjoint(nullspacemat(tra2)) * x
Δq = UnitQuaternion(rand(4)...)
setPosition!(body1, body2; p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δx = Δx, Δq = Δq)
@test norm(minimalCoordinates(tra2, body1, body2) - x[1], Inf) < 1e-10

v = srand(1)
Δv = zerodimstaticadjoint(nullspacemat(tra2)) * v
Δω = rand(3)
setVelocity!(body1, body2; p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δv = Δv, Δω = Δω)
@test norm(minimalVelocities(tra2, body1, body2) - v[1], Inf) < 1e-10


################################################################################
# Test min -> max -> min
################################################################################

# hopper
mech = getmechanism(:hopper)
Random.seed!(100)
x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(2)]
# x0 = [zeros(3); vector(UnitQuaternion(1,0.1,0,0)); zeros(3); 0.0; 0.2; 0.0; zeros(2)]
z0 = min2max(mech, x0)
x1 = max2min(mech, z0)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
@test norm(x0[14:15] - x1[14:15], Inf) < 1e-10


g0 = -9.81
Nlink0 = 2
mech = getmechanism(:snake, g = g0, Nlink = Nlink0, jointtype = :Spherical, damper = 1.0)
initialize!(mech, :snake, ω = [0.1, 0.2, 0.0])

function controller!(mechanism, k)
	# for (i,eqc) in enumerate(collect(mechanism.eqconstraints))
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints[2:end]))
        nu = controldim(eqc)
        setForce!(mechanism, eqc, 0.1*SVector{nu}(Δt0 * (rand(nu) .- 2.0) ))
    end
    return
end
@elapsed storage = simulate!(mech, 3.00, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = 1e-6)
visualize(mech, storage, vis = vis)




Random.seed!(100)
nx = minCoordDim(mech)
x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); zeros(nx - 13)]
# x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(3); 0.0; 0.0; 0.1]
# x0 = [zeros(3); vector(UnitQuaternion(1,0,0,0.)); zeros(3); zeros(3); rand(3); 0.0; 0.0; 0.1]
# x0 = [zeros(3); vector(UnitQuaternion(1,0,0,0.)); zeros(3); zeros(3); 0.1; 0.0; 0.0; 0.0; 0.1; 0.0]
z0 = min2max(mech, x0)
storage = Storage(2,2)
for t = 1:2
	for b = 1:2
		x2, v15, q2, ϕ15 = unpackMaxState(z0, b)
		storage.x[b][t] = x2
		storage.v[b][t] = v15
		storage.q[b][t] = q2
		storage.ω[b][t] = ϕ15
	end
end

visualize(mech, storage, vis = vis)

x1 = max2min(mech, z0)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
@test norm(x0[14:16] - x1[14:16], Inf) < 1e-10
@test norm(x0[17:nx] - x1[17:nx], Inf) < 1e-10


# snake
for jointtype in jointtypes
	@show jointtype

	Nlink0 = 2
	mech = getmechanism(:snake, Nlink = Nlink0, jointtype = jointtype)
	Random.seed!(100)
	nx = minCoordDim(mech)
	x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx - 13)]
	z0 = min2max(mech, x0)
	x1 = max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
	@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10
end

for jointtype in jointtypes
	@show jointtype

	Nlink0 = 2
	mech = getmechanism(:twister, Nlink = Nlink0, jointtype = jointtype)
	Random.seed!(100)
	nx = minCoordDim(mech)
	x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx - 13)]
	z0 = min2max(mech, x0)
	x1 = max2min(mech, z0)
	@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
	@test quaterror(x0[4:7], x1[4:7]) < 1e-10
	@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
	@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
	@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10
end

# humanoid
mech = getmechanism(:humanoid)
Random.seed!(100)
nx = minCoordDim(mech)
x0 = [rand(3); vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx - 13)]
z0 = min2max(mech, x0)
x1 = max2min(mech, z0)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10


vis = Visualizer()
open(vis)







# atlas
mech = getmechanism(:atlas, model_type = :simple, contact = true, damper = 10.0)
# mech = getmechanism(:atlas, model_type = :debug, contact = false)


# body1 = collect(mech.bodies)[1]
# body2 = collect(mech.bodies)[2]
# body3 = collect(mech.bodies)[3]
#
# eqc1 = collect(mech.eqconstraints)[1]
# eqc2 = collect(mech.eqconstraints)[2]
# eqc3 = collect(mech.eqconstraints)[3]
#
# pbody1 = getbody(mech, eqc1.parentid)
# pbody2 = getbody(mech, eqc2.parentid)
# pbody3 = getbody(mech, eqc3.parentid)
#
# cbody1 = getbody(mech, eqc1.childids[1])
# cbody2 = getbody(mech, eqc2.childids[1])
# cbody3 = getbody(mech, eqc3.childids[1])
#
#
# mech.origin# origin
# eqc1
# body2# pelvis
# eqc3
# body1# ltorso
# eqc2
# body3# mtorso

# zinit = getMaxState(mech)
# xinit = max2min(mech, zinit)
# xinit[1:13]
# norm(xinit[14:end], Inf)
# vizMax(mech, zinit)

Random.seed!(10)
nx = minCoordDim(mech)
x0 = [[0,0,1+rand()]; vector(UnitQuaternion(rand(4)...)); rand(3); rand(3); rand(nx-13)]
# x0 = [zeros(3); vector(UnitQuaternion(1,0,0,0)); zeros(3); zeros(3); ones(nx-13)]
# x0 = [zeros(3); vector(UnitQuaternion(1,0,0,0)); zeros(3); zeros(3); ones(2); ones(2)]
# x0 = [zeros(3); vector(UnitQuaternion(1,0,0,0)); zeros(3); zeros(3); 1.0rand(nx-13)]

z0 = min2max(mech, deepcopy(x0))

@elapsed storage = simulate!(mech, 8.0, controller!, record = true,
	solver = :mehrotra!, verbose = false, ϵ = 1e-4)
visualize(mech, storage, vis = vis)



# norm(x0 - xinit, Inf)
# norm(z0 - zinit, Inf)

vizMax(deepcopy(mech), deepcopy(z0))
x1 = max2min(deepcopy(mech), deepcopy(z0))
# z1 = min2max(mech, x1)
@test norm(x0[1:3] - x1[1:3], Inf) < 1e-10
@test quaterror(x0[4:7], x1[4:7]) < 1e-10
@test norm(x0[8:10] - x1[8:10], Inf) < 1e-10
@test norm(x0[11:13] - x1[11:13], Inf) < 1e-10
@test norm(x0[14:nx] - x1[14:nx], Inf) < 1e-10

plot(x0[14:nx] - x1[14:nx])



# pendulum
mech = getmechanism(:pendulum)
Random.seed!(100)
x0 = rand(2)
z0 = min2max(mech, x0)
x1 = max2min(mech, z0)
@test norm(x0 - x1, Inf) < 1e-10

# Npendulum
for jointtype in jointtypes
	@show jointtype
	Nlink0 = 2
	mech = getmechanism(:npendulum, Nlink = Nlink0, jointtype = jointtype)
	Random.seed!(100)
	nx = minCoordDim(mech)
	x0 = rand(nx)
	z0 = min2max(mech, x0)
	x1 = max2min(mech, z0)
	@test norm(x0 - x1, Inf) < 1e-10
end


# ################################################################################
# # Test max -> min -> max
# ################################################################################
#
# vis = Visualizer()
# open(vis)
#
# Δt0 = 0.01
# g0 = -0.00
# spring0 = 0.0
# damper0 = 1.0
# mech = getmechanism(:hopper, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
# initialize!(mech, :hopper, v = [0.1, 0.2, 0.3], ω = [0.4, 0.4, 0.1])
#
# @elapsed storage = simulate!(mech, 3.00, record = true, solver = :mehrotra!, verbose = false, ϵ = 1e-12)
# visualize(mech, storage, vis = vis)
#
# ################################################################################
# # Test setPosition! and setVelocity!
# ################################################################################
# Nb = length(storage.x)
# Nt = length(storage.x[1])
# err = []
# for t = 1:Nt-1
# 	z0 = zeros(13Nb)
# 	for i = 1:Nb
# 		x2 = storage.x[i][t]
# 		# v15 = storage.v[i][t] * 1.000000000000000000
# 		v15 = storage.vl[i][t] * 1.000000000000000000
# 		q2 = storage.q[i][t]
# 		# ϕ15 = storage.ω[i][t] * 1.0000000000000000000
# 		ϕ15 = storage.ωl[i][t] * 1.0000000000000000000
# 		z0[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
# 	end
# 	x0 = max2min(mech, z0)
# 	z1 = min2max(mech, x0)
#
# 	# push!(err, norm((z0 - z1)[1:3], Inf))
# 	# push!(err, norm((z0 - z1)[4:6], Inf))
# 	# push!(err, norm(quaterror(z0[7:10], z1[7:10])))
# 	# push!(err, norm((z0 - z1)[11:13], Inf))
# 	# push!(err, norm((z0 - z1)[13 .+ (1:3)], Inf))
# 	push!(err, norm((z0 - z1)[13 .+ (4:6)], Inf))
# 	# push!(err, norm(quaterror(z0[13 .+ (7:10)], z1[13 .+ (7:10)])))
# 	# push!(err, norm((z0 - z1)[13 .+ (11:13)], Inf))
# end
# norm(err, Inf)
# plot(log.(10,err))
