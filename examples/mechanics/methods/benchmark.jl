function benchmark_momentum(timestep::Vector; 
    Nsim=1, 
    tsim=1.0, 
    tctrl=1.0, 
    seed=0, 
    ϵ=1.0e-14, 
    control_amplitude=0.05)

    Nt = length(timestep)
    speed = zeros(Nt)
    plin = zeros(Nt)
    pang = zeros(Nt)
    for i = 1:Nt
        @show timestep[i]
        mech, storage, tcompute = astronaut_simulation(
            Nsim=Nsim, 
            timestep=timestep[i],
            tsim=tsim, 
            tctrl=tctrl, 
            seed=seed, 
            ϵ=ϵ, 
            control_amplitude=control_amplitude)
        speed[i], plin[i], pang[i] = process_momentum_runs(mech, storage, tcompute)
    end
    return speed, plin, pang
end

function benchmark_energy(timestep::Vector; 
    Nsim=1, 
    tsim=2.0, 
    tctrl=1.0, 
    seed=0, 
    ϵ=1.0e-14, 
    control_amplitude=0.05)

    Nt = length(timestep)
    speed = zeros(Nt)
    energy = zeros(Nt)
    for i = 1:Nt
        @show timestep[i]
        mech, storage, tcompute = astronaut_simulation(
            Nsim=Nsim, 
            timestep=timestep[i],
            tsim=tsim, 
            tctrl=tctrl, 
            seed=seed, 
            ϵ=ϵ, 
            control_amplitude=control_amplitude)
        speed[i], energy[i] = process_energy_runs(mech, storage, tcompute, tctrl)
    end
    return speed, energy
end