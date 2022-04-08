using BenchmarkTools
using Plots

################################################################################
# Query NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()
point = rand(Float32, 100,3) ./ 5
density = py"density_query"(nerf_object, point)
density_grad = py"density_gradient_query"(nerf_object, point)
@belapsed density = py"density_query"(nerf_object, point)
@belapsed density_grad = py"density_gradient_query"(nerf_object, point)

# evaluation & gradient benchmark
Ns = [1,3,10,30,100,300,1000,3000,10000]
points = [rand(Float32, N,3) for N in Ns]
eval_times = [@elapsed py"density_query"(nerf_object, point) for point in points]
grad_times = [@elapsed py"density_gradient_query"(nerf_object, point) for point in points]
plot(Ns, eval_times, axis=:log, label="evaluation", xlabel="number of samples", ylabel="time (s)",
    title="batch querying is more efficient")
plot!(Ns, grad_times, axis=:log, label="gradient", xlabel="number of samples", ylabel="time (s)",
    title="batch querying is more efficient")
