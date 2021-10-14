include("../test/example_file_names.jl")

for i=1:length(files)
    include("../test/examples/"*files[i]*".jl")
    mech.g = 0.0

    steps = Base.OneTo(100)
    storage = Storage{Float64}(steps,length(mech.bodies))

    if files[i] ∈ controlled
        control_function = eval(Meta.parse(files[i]*"_control!"))
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $control_function) samples=2
    else
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage) samples=2
    end
end
