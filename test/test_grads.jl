###############################################################################
# Particle
################################################################################

mech = getmechanism(:box, mode = :box, contact = true)
initialize!(mech, :box, x = [0,0,0.], q = one(UnitQuaternion), v = [0,0,0.], ω = [0,0,0.])
initialize!(mech, :box, x = rand(3), q = UnitQuaternion(1,0,0,0.), v = rand(3), ω = rand(3))
initialize!(mech, :box, x = 1 .+ rand(3), q = UnitQuaternion(rand(4)...), v = rand(3), ω = rand(3))

nx = minCoordDim(mech)
nz = maxCoordDim(mech)
nu = controldim(mech)

z0 = getMaxState(mech)
x0 = max2min(mech, z0)
u0 = mech.Δt * rand(nu)

z1 = step!(mech, z0, u0, ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)
x1 = max2min(mech, z1)

∇x0 = ∇min2max(mech, x0) # dz(x)/dx
∇z0 = ∇max2min(mech, z0) # dx(z)/dz
@test norm(∇z0 * ∇x0 - I(nx)) < 1e-7

∇z_z, ∇u_z = getMaxGradients!(mech, z0, u0)
∇x_x, ∇u_x = getMinGradients!(mech, z0, u0)

FDz_z = FiniteDiff.finite_difference_jacobian(z -> step!(mech, z, u0,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false), z0)
FDx_x = FiniteDiff.finite_difference_jacobian(x -> max2min(mech, step!(mech, min2max(mech,x), u0,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)), x0)

FDu_z = FiniteDiff.finite_difference_jacobian(u -> step!(mech, z0, u,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false), u0)
FDu_x = FiniteDiff.finite_difference_jacobian(u -> max2min(mech, step!(mech, min2max(mech,x0), u,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)), u0)

norm(FDz_z - ∇z_z, Inf)
norm(FDx_x - ∇x_x, Inf)
plot(Gray.(FDx_x))
plot(Gray.(∇x_x))
plot(Gray.(abs.(1e4(FDx_x - ∇x_x))))

norm(FDu_z - ∇u_z, Inf)
norm(FDu_x - ∇u_x, Inf)
plot(Gray.(FDu_x))
plot(Gray.(∇u_x))
plot(Gray.(abs.(1e4(FDu_x - ∇u_x))))





################################################################################
# Hopper
################################################################################

mech = getmechanism(:hopper, contact = false)
initialize!(mech, :hopper, altitude = 1.0)

nx = minCoordDim(mech)
nz = maxCoordDim(mech)
nu = controldim(mech)

z0 = getMaxState(mech)
x0 = max2min(mech, z0)
u0 = mech.Δt * rand(nu)

z1 = step!(mech, z0, u0, ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)
x1 = max2min(mech, z1)

∇x0 = ∇min2max(mech, x0) # dz(x)/dx
∇z0 = ∇max2min(mech, z0) # dx(z)/dz
@test norm(∇z0 * ∇x0 - I(nx)) < 1e-7

	#max #min
norm(∇x0[1:3, 1:3] - I(3), Inf) < 1e-8
norm(∇x0[4:6, 8:10] - I(3), Inf) < 1e-8
norm(∇x0[7:10, 4:7] - Diagonal([0,1,1,1]), Inf) < 1e-8
norm(∇x0[11:13, 11:13] - I(3), Inf) < 1e-8

∇z_z, ∇u_z = getMaxGradients!(mech, z0, u0)
∇x_x, ∇u_x = getMinGradients!(mech, z0, u0)

FDz_z = FiniteDiff.finite_difference_jacobian(z -> step!(mech, z, u0,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false), z0)
FDx_x = FiniteDiff.finite_difference_jacobian(x -> max2min(mech, step!(mech, min2max(mech,x), u0,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)), x0)

FDu_z = FiniteDiff.finite_difference_jacobian(u -> step!(mech, z0, u,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false), u0)
FDu_x = FiniteDiff.finite_difference_jacobian(u -> max2min(mech, step!(mech, min2max(mech,x0), u,
	ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)), u0)

norm(FDz_z - ∇z_z, Inf)
norm(FDx_x - ∇x_x, Inf)
plot(Gray.(FDx_x))
plot(Gray.(∇x_x))
plot(Gray.(abs.(FDx_x - ∇x_x)))
FDx_x[11:13,4:7]
∇x_x[11:13,4:7]
(FDx_x - ∇x_x)[11:13,4:7]
FDx_x[10:12,4:6]



∇x_x[10:12,4:6]



(FDx_x - ∇x_x)[10:12,4:6]



plot(Gray.(FDz_z))
plot(Gray.(∇z_z))
plot(Gray.(abs.(FDz_z - ∇z_z)))

norm(FDu_z - ∇u_z, Inf)
norm(FDu_x - ∇u_x, Inf)
plot(Gray.(FDu_x))
plot(Gray.(∇u_x))
plot(Gray.(abs.(FDu_x - ∇u_x)))
