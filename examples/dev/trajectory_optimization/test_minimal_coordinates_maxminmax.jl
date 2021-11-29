# vis = Visualizer()
# open(vis)

Δt0 = 0.01
g0 = -0.00
spring0 = 0.0
damper0 = 1.0
mech = getmechanism(:hopper, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :hopper, v = [0.1, 0.2, 0.3], ω = [0.4, 0.4, 0.1])

@elapsed storage = simulate!(mech, 3.00, record = true, solver = :mehrotra!, verbose = false, ϵ = 1e-12)
visualize(mech, storage, vis = vis)

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
# Test setPosition! and setVelocity!
################################################################################
Nb = length(storage.x)
Nt = length(storage.x[1])
err = []
for t = 1:Nt-1
	z0 = zeros(13Nb)
	for i = 1:Nb
		x2 = storage.x[i][t]
		# v15 = storage.v[i][t] * 1.000000000000000000
		v15 = storage.vl[i][t] * 1.000000000000000000
		q2 = storage.q[i][t]
		# ϕ15 = storage.ω[i][t] * 1.0000000000000000000
		ϕ15 = storage.ωl[i][t] * 1.0000000000000000000
		z0[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
	end
	x0 = max2min(mech, z0)
	z1 = min2max(mech, x0)

	# push!(err, norm((z0 - z1)[1:3], Inf))
	# push!(err, norm((z0 - z1)[4:6], Inf))
	# push!(err, norm(quaterror(z0[7:10], z1[7:10])))
	# push!(err, norm((z0 - z1)[11:13], Inf))
	# push!(err, norm((z0 - z1)[13 .+ (1:3)], Inf))
	push!(err, norm((z0 - z1)[13 .+ (4:6)], Inf))
	# push!(err, norm(quaterror(z0[13 .+ (7:10)], z1[13 .+ (7:10)])))
	# push!(err, norm((z0 - z1)[13 .+ (11:13)], Inf))
end
norm(err, Inf)
plot(log.(10,err))
plot([norm(v) for v in storage.v[2]])
plot(err[3744:3750])

t0 = 3744
z0 = zeros(13Nb)
for i = 1:Nb
	x2 = storage.x[i][t0] * 1.0
	# v15 = storage.v[i][t0] * 1.0
	v15 = storage.vl[i][t0] * 1.0
	q2 = storage.q[i][t0] * 1.0
	# ϕ15 = storage.ω[i][t0] * 1.0
	ϕ15 = storage.ωl[i][t0] * 1.0
	z0[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
end

z0[13 .+ (1:3)] -= z0[1:3]
z0[1:3] = zeros(3)
z0[13 .+ (4:6)] -= z0[4:6]
z0[4:6] = zeros(3)
q1 = UnitQuaternion(z0[7:10]...)
q2 = UnitQuaternion(z0[13 .+ (7:10)]...)

z0[7:10] = vector(one(UnitQuaternion))
z0[13 .+ (7:10)] = vector(inv(q1) * q2)
z0[1:3] = vrotate(z0[1:3], inv(q1))
z0[4:6] = vrotate(z0[4:6], inv(q1))
z0[11:13] = z0[11:13]
z0[13 .+ (1:3)] = vrotate(z0[13 .+ (1:3)], inv(q1))
z0[13 .+ (4:6)] = vrotate(z0[13 .+ (4:6)], inv(q1))

q1n = UnitQuaternion(z0[7:10]...)
q2n = UnitQuaternion(z0[13 .+ (7:10)]...)

ω1w = vrotate(z0[11:13], q1n)
ω2w = vrotate(z0[13 .+ (11:13)], q2n)
ω1w - ω2w

z0[13 .+ (4:6)] -= skew(ω1w) * (z0[13 .+ (1:3)] - z0[1:3])
z0[11:13] = zeros(3)
z0[13 .+ (11:13)] -= vrotate(ω1w, inv(q2n))
ω1w = vrotate(z0[11:13], q1n)
ω2w = vrotate(z0[13 .+ (11:13)], q2n)
ω1w - ω2w

@show z0

x0 = max2min(mech, z0)
z1 = min2max(mech, x0)
norm((z0 - z1)[(4:6)], Inf)
norm((z0 - z1)[13 .+ (4:6)])
z0[13 .+ (4:6)]



z1[13 .+ (4:6)]



@show z0[1:13]
@show z0[13 .+ (1:13)]



# Max coord by hand: translational
x21 = zeros(3)
q21 = UnitQuaternion(1,0,0,0.)
v151 = zeros(3)
ϕ151 = [0.1, 0.2, 0.0]

x22 = [0,0,-0.5]
q22 = UnitQuaternion(1,0,0,0.)
v152 = [-0.1,0.05,-0.2]
ϕ152 = [0.1, 0.2, 0.0]

z0 = [x21; v151; vector(q21); ϕ151;
	  x22; v152; vector(q22); ϕ152]
x0 = max2min(mech, z0)
norm(x0[1:3] - z0[1:3], Inf)
norm(x0[8:10] - z0[4:6], Inf)
quaterror(x0[4:7], z0[7:10])
norm(x0[11:13] - z0[11:13], Inf)
norm(x0[14] + 0.5, Inf)
norm(x0[15] + 0.2, Inf)

z1 = min2max(mech, x0)
norm(z0 - z1, Inf)
norm((z0 - z1)[1:3], Inf)
norm((z0 - z1)[4:6], Inf)
norm(quaterror(z0[7:10], z1[7:10]))
norm((z0 - z1)[11:13], Inf)
norm((z0 - z1)[13 .+ (1:3)], Inf)
norm((z0 - z1)[13 .+ (4:6)], Inf)
norm(quaterror(z0[13 .+ (7:10)], z1[13 .+ (7:10)]))
norm((z0 - z1)[13 .+ (11:13)], Inf)


# Max coord by hand: translational
x21 = zeros(3)
q21 = UnitQuaternion(RotX(π))
v151 = zeros(3)
ϕ151 = [0.1, 0.2, 0.0]
# ϕ151 = [0.1, 0.0, 0.0]

x22 = [0,0,0.5]
q22 = UnitQuaternion(RotX(π))
v152 = [-0.1,-0.05,0.2]
# v152 = [0,-0.05,0.2]
ϕ152 = [0.1, 0.2, 0.0]
# ϕ152 = [0.1, 0.0, 0.0]

# Max coord by hand: translational
x21 = zeros(3)
q21 = UnitQuaternion(RotX(π/2))
v151 = zeros(3)
ϕ151 = [0.1, 0.2, 0.0]
# ϕ151 = [0.1, 0.0, 0.0]

x22 = [0,0.5,0]
q22 = UnitQuaternion(RotX(π/2))
v152 = [-0.1,0.2,0.05]
# v152 = [0,0.2,0.05]
ϕ152 = [0.1, 0.2, 0.0]
# ϕ152 = [0.1, 0.0, 0.0]

# Max coord by hand: translational
x21 = zeros(3)
q21 = UnitQuaternion(RotX(π/4))
v151 = zeros(3)
# ϕ151 = [0.1, 0.2, 0.0]
ϕ151 = [0.1, 0.0, 0.0]

x22 = [0,-sqrt(2)/4,sqrt(2)/4]
q22 = UnitQuaternion(RotX(π/4))
# v152 = [-0.1,0.2,0.05]
v152 = [0, -0.05*sqrt(2)/2, -0.05*sqrt(2)/2] + [0, -0.2*sqrt(2)/2, 0.2*sqrt(2)/2]
# ϕ152 = [0.1, 0.2, 0.0]
ϕ152 = [0.1, 0.0, 0.0]


z0 = [x21; v151; vector(q21); ϕ151;
	  x22; v152; vector(q22); ϕ152]
x0 = max2min(mech, z0)
norm(x0[1:3] - z0[1:3], Inf)
norm(x0[8:10] - z0[4:6], Inf)
quaterror(x0[4:7], z0[7:10])
norm(x0[11:13] - z0[11:13], Inf)
norm(x0[14] - 0.5, Inf)
norm(x0[15] - 0.2, Inf)
x0[14]
x0[15]

z1 = min2max(mech, x0)
norm(z0 - z1, Inf)
norm((z0 - z1)[1:3], Inf)
norm((z0 - z1)[4:6], Inf)
norm(quaterror(z0[7:10], z1[7:10]))
norm((z0 - z1)[11:13], Inf)
norm((z0 - z1)[13 .+ (1:3)], Inf)
norm((z0 - z1)[13 .+ (4:6)], Inf)
norm(quaterror(z0[13 .+ (7:10)], z1[13 .+ (7:10)]))
norm((z0 - z1)[13 .+ (11:13)], Inf)

z0[13 .+ (4:6)]
z1[13 .+ (4:6)]

z0[13 .+ (1:3)]
z1[13 .+ (1:3)]
