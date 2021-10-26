using BenchmarkTools

include("example_file_names.jl")

for file in files
    include("examples/"*file*".jl")
    if file ∈ controlled
        storage = simulate!(mech, 10., eval(Meta.parse(file*"_control!")), record = true)
    else
        storage = simulate!(mech, 10., record = true, solver = :mehrotra!, verbose = false)
    end
    n = length(storage.x)
    for i=1:n
        @test !any(ConstrainedDynamics.sisnan.(storage.x[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.q[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.v[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.ω[i]))
    end
end

include("examples/dummycontroller.jl")
storage = simulate!(mech, 10., dummycontroller_controller, record = true, debug = true)
n = length(storage.x)
for i=1:n
    @test !any(ConstrainedDynamics.sisnan.(storage.x[i]))
    @test !any(ConstrainedDynamics.sisnan.(storage.q[i]))
    @test !any(ConstrainedDynamics.sisnan.(storage.v[i]))
    @test !any(ConstrainedDynamics.sisnan.(storage.ω[i]))
end
