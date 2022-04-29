# include("../../robodojo/halfcheetah/model.jl")
# include("../../robodojo/halfcheetah/visuals.jl")
# include("../../robodojo/dynamics.jl")
#
# RoboDojo.RESIDUAL_EXPR
# force_codegen = true
# # force_codegen = false
# robot = halfhyena
# include("../../robodojo/codegen.jl")
# RoboDojo.RESIDUAL_EXPR








# du0 = zeros(nx + nθ, nx + nθ)
# x0 = ones(nx)
# u0 = ones(nx + nθ)
# w0 = ones(nw)
# f1u(du0, x0, u0, w0)
# @benchmark $f1u($du0, $x0, $u0, $w0)
#
# dx0 = zeros(nx + nθ, nx + nθ)
# x0 = ones(nx + nθ)
# u0 = ones(nx)
# w0 = ones(nw)
# ftx(dx0, x0, u0, w0)
# @benchmark $ftx($dx0, $x0, $u0, $w0)


# dx[1:nq,1:nq]
# dx0[1:nq,1:nq]
# dx[1:nq,1:nq] - dx0[1:nq,1:nq]
# norm(dx[1:nq,1:nq] - dx0[1:nq,1:nq])
#
# dx[1:nq,nq .+ (1:nq)]
# dx0[1:nq,nq .+ (1:nq)]
# dx[1:nq,nq .+ (1:nq)] - dx0[1:nq,nq .+ (1:nq)]
# norm(dx[1:nq,nq .+ (1:nq)] - dx0[1:nq,nq .+ (1:nq)])
#
# dx[nq .+ (1:nq),1:nq]
# dx0[nq .+ (1:nq),1:nq]
# dx[nq .+ (1:nq),1:nq] - dx0[nq .+ (1:nq),1:nq]
# norm(dx[nq .+ (1:nq),1:nq] - dx0[nq .+ (1:nq),1:nq])
#
#
#
# plot(Gray.(dx))
# plot(Gray.(dx0))
# plot(Gray.(dx - dx0))
# plot(Gray.(1e2abs.(dx - dx0)))
# plot(Gray.(1e4abs.(dx - dx0)))
