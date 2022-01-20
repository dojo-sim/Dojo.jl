@testset "Behaviors: Quadruped" begin
    mech = getmechanism(:quadruped, Δt=0.05, g=-9.81, cf=0.8, damper=1000.0, spring=30.0)
    initialize!(mech, :quadruped)
    try
        storage = simulate!(mech, 5.0, record=true, verbose=false)
        @test true
    catch 
        @test false 
    end
end

@testset "Behavior: Box toss" begin
    for Δt in [0.10, 0.05, 0.01, 0.005]
        mech = getmechanism(:box, Δt = Δt, g = -9.81, cf = 0.1)
        initialize!(mech, :box, x=[0,0,0.5], v=[1,1.5,1.], ω=[5,4,2.] .* Δt)
        storage = simulate!(mech, 5.0, record=true,
            opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
        # visualize(mech, storage, vis=vis)
        @test norm(storage.v[1][end], Inf) < 1e-12
        @test norm(storage.x[1][end][3] - 0.25, Inf) < 1e-3
    end
end
