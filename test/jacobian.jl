function test_solmat(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, ctrl::Any=(m,k)->nothing,
        Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "solmat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g; kwargs...)
        initialize!(mechanism, model)

        storage = simulate!(mechanism, tsim, ctrl, record = true, solver = :mehrotra!, verbose = false)
        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        solmat = full_matrix(mechanism.system)
        # finite diff
        fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ = 1e-5, verbose = verbose)
        @test norm(fd_solmat + solmat, Inf) < ϵ
    end
    return nothing
end

function test_datamat(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, ctrl::Any=(m,k)->nothing,
        Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "datamat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g; kwargs...)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, tsim, ctrl, record = true, solver = :mehrotra!, verbose = false)

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        datamat0 = full_data_matrix(mechanism, attjac = true)
        datamat1 = full_data_matrix(mechanism, attjac = false)
        # finite diff
        fd_datamat1 = finitediff_data_matrix(mechanism, data, sol, δ = 1e-5, verbose = verbose)
        fd_datamat0 = fd_datamat1 * attjac

        @test norm((fd_datamat0 + datamat0), Inf) < ϵ
        @test norm((fd_datamat1 + datamat1), Inf) < ϵ
    end
    return nothing
end

function test_sensitivity(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, ctrl::Any=(m,k)->nothing,
        Δt::T = 0.01, g::T = -9.81, cf::T = 0.8,
        contact::Bool = true, verbose::Bool = false) where {T}

    @testset "sensitivity: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g, cf = cf, contact = contact)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, tsim, ctrl, record = true, solver = :mehrotra!, verbose = false)

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        datamat = full_data_matrix(mechanism)
        solmat = full_matrix(mechanism.system)
        sensi = - (solmat \ datamat)

        # finite diff
        fd_sensi = finitediff_sensitivity(mechanism, data, verbose = verbose) * attjac
        @test norm(fd_sensi - sensi) / norm(fd_sensi) < ϵ
    end
    return nothing
end

function cont!(mechanism, k; u = 0.1)
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        nu = controldim(eqc, ignore_floating_base = false)
        su = mechanism.Δt * u * sones(nu)
        setForce!(mechanism, eqc, su)
    end
    return
end

#TODO: decrease ϵ tolerance once we replace finite-difference methods
# Flying after 0.1 sec simulation
tsim = 0.1
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,      tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)

tsim = 0.1
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,      tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)

# In contact with the ground after 0.4 sec simulation
tsim = 0.4
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,      tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)

tsim = 0.4
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,      tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)

function cont!(mechanism, k; u = 0.1)
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        nu = controldim(eqc, ignore_floating_base = false)
        su = mechanism.Δt * u * sones(nu)
        setForce!(mechanism, eqc, su)
    end
    return
end

function plot_cone(s, Δs; plt = plot())
	θ = 0:0.03:2π
	Xc = cos.(θ)
	Yc = sin.(θ)
	plot!(Xc, Yc, legend = false, aspectratio = 1.0, xlims = (-2,2), ylims = (-2,2))
	scatter!(1/s[1]*s[2:2], 1/s[1]*s[3:3], markersize = 6.0)
	α = 0:0.05:1.0
	Sα = [s + αi*Δs for αi in α]
	Xα = [1/s[1]*s[2] for s in Sα]
	Yα = [1/s[1]*s[3] for s in Sα]
	plot!(Xα, Yα)
	display(plt)
	return nothing
end


# s  = [1.001619650763305e+0, 1.737348210974044e-1, -1.316707656377632e-1]
# Δs = [9.525076266034363e+2, 8.050194800248429e+0,  6.473898126907718e+0]
# s  = [1.001619650763305e+0, 1.737348210974044e-1, -1.316707656377632e-1]
# Δs = [8.945577700909580e+2, 7.678512219876204e+0,  5.985574831920998e+0]
# sp  =[ 1.001619650763305e+0,  1.737348210974044e-1, -1.316707656377632e-1]
# Δsp =[ 9.525076266034363e+2,  8.050194800248429e+0,  6.473898126907718e+0]
# sc  =[ 1.001619650763305e+0,  1.737348210974044e-1, -1.316707656377632e-1]
# Δsc =[ 8.945577700909580e+2,  7.678512219876204e+0,  5.985574831920998e+0]

# γp  = [2.0806866959138824e-15, -3.5160872884030105e-7, -4.196336237276323e-7]
# Δγp = [-2.93464026390496e-15, 0.0003347200444208716, 0.00039947752616968673]
# γc  = [2.0806866959138824e-15, -3.5160872884030105e-7, -4.196336237276323e-7]
# Δγc = [-1.681677547721007e-15, -0.3179934702225138, -0.3795149019007642]
# plt = plot()
# plot_cone(γp, Δγp, plt = plt)
# plot_cone(γc, Δγc, plt = plt)
