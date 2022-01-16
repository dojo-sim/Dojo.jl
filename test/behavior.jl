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
