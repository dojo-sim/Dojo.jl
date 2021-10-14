using Test

@testset "Simulation" begin
    # Atlas
    mech = getmechanism(:atlas, contact = false)
    initialize!(mech, :atlas)
    storage = simulate!(mech, 0.2, record = true, solver = :mehrotra!)

    # Dice Linear
    linmech = getmechanism(:dice, conetype = :linear)
    initialize!(linmech, :dice)
    linstorage = simulate!(linmech, 0.5, record = true, solver = :mehrotra!)
    # Dice SOC
    socmech = getmechanism(:dice, conetype = :soc)
    initialize!(socmech, :dice)
    socstorage = simulate!(socmech, 0.5, record = true, solver = :mehrotra!)

    # Snake Linear
    linmech = getmechanism(:snake, conetype = :linear)
    initialize!(linmech, :snake)
    linstorage = simulate!(linmech, 0.5, record = true, solver = :mehrotra!)
    # Snake SOC
    socmech = getmechanism(:snake, conetype = :soc)
    initialize!(socmech, :snake)
    socstorage = simulate!(socmech, 0.5, record = true, solver = :mehrotra!)

    # Npendulum
    mech = getmechanism(:npendulum)
    initialize!(mech, :npendulum)
    storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)

    # Pendulum
    mech = getmechanism(:pendulum)
    initialize!(mech, :pendulum)
    storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)
end
