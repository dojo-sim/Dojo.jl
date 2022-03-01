using Symbolics

nx = 6
nu = 1
@variables δx[1:nx] δu[1:nu] Bf[1:nx*nu]

function cost(δx, δu, Bf)
    B = reshape(Bf, nx, nu)
    return 1/2 * (δx - B*δu)' * (δx - B*δu)
end

c = cost(δx, δu, Bf)
dc = Symbolics.gradient(c, Bf)
ddc = Symbolics.hessian(c, Bf)

c_fct = eval(Symbolics.build_function([c], δx, δu, Bf)[1])
dc_fct = eval(Symbolics.build_function(dc, δx, δu, Bf)[1])
ddc_fct = eval(Symbolics.build_function(ddc, δx, δu, Bf)[1])

δx0 = rand(nx)
δu0 = rand(nu)
Bf0 = rand(nx*nu)

c_fct(δx0, δu0, Bf0)
dc_fct(δx0, δu0, Bf0)
ddc_fct(δx0, δu0, Bf0)

function leastsquares(δx, δu)
    N = length(δx)
    nx = length(δx[1])
    nu = length(δu[1])
    Bf = zeros(nx*nu)
    e = 0.0
    grad = zeros(nx*nu)
    hess = zeros(nx*nu, nx*nu)
    for i = 1:N
        e += c_fct(δx[i], δu[i], Bf)[1]
        grad += dc_fct(δx[i], δu[i], Bf)
        hess += ddc_fct(δx[i], δu[i], Bf)
    end
    Bf = Bf - hess \ grad
    B = reshape(Bf, (nx, nu))
    return B
end