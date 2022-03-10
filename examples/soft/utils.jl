# ENV["PYCALL_JL_RUNTIME_PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"
# ENV["PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"

using BenchmarkTools
using Plots
using Pkg
using PyCall
# Pkg.build("PyCall")

osf_path = joinpath("/home/simon/research/repos/osf-pytorch")
# pushfirst!(pyimport("sys")."path", "")
pushfirst!(pyimport("sys")."path", osf_path)
@pyinclude(joinpath(osf_path, "extract_density_julia.py"))


################################################################################
# Query NeRf density and its gradient
################################################################################
function generate_test_nerf()
    render_kwargs_test = py"generate_test_nerf"()
    return render_kwargs_test
end

function density_query(render_kwargs_test::Dict, p::AbstractMatrix)
    density = py"density_query"(render_kwargs_test, p)
    return density
end

function density_gradient_query(render_kwargs_test::Dict, p::AbstractMatrix)
    density_grad = py"density_gradient_query"(render_kwargs_test, p)
    return density_grad
end


# p = rand(Float32, 100,3) ./ 5
# render_kwargs_test = generate_test_nerf()
# density = density_query(render_kwargs_test, p)
# density_grad = density_gradient_query(render_kwargs_test, p)
