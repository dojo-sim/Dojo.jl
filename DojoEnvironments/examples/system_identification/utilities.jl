# Taken from Nocedal and Wright, Algorithm 6.1
function quasi_newton_solve(f, fgH, x0; ftol=-Inf, gtol=1e-4, iter=100, α0=1.0,
        lower_bound=-Inf, upper_bound=Inf, reg=1e-3, reg_min=1e-9, reg_max=1e6)

    x = copy(x0)
    X = [copy(x0)]
    ls_failure = false

    for k = 1:iter
        fe, ge, He = fgH(x)
        He += reg * I
        if ls_failure
            reg = clamp(reg * 2, reg_min, reg_max)
        else
            reg = clamp(reg / 1.5, reg_min, reg_max)
        end
        ((norm(ge, Inf) < gtol) || (fe < ftol)) && break
        p = - He \ ge
        α, ls_failure = clamped_linesearch(f, x, p, fe; α0=α0, lower_bound, upper_bound)
        x = clamp.(x + α*p, lower_bound, upper_bound)
        push!(X, copy(x))
        println("k:", k, "   f:", Dojo.scn(fe, digits=3))
    end
    return x, X
end

function clamped_linesearch(f, x, p, fprev; α0=1.0, iter=4,
        lower_bound=-Inf, upper_bound=Inf)
    α = α0
    ls_failure = false
    for k = 1:iter
        xc = clamp.(x + α*p, lower_bound, upper_bound)
        (f(xc) <= fprev) && break
        α /= 3
        (k == iter) && (α = 0.001 / norm(p,Inf); ls_failure = true)
    end
    return α, ls_failure
end



function get_contact_gradients!(mechanism::Mechanism, z::AbstractVector, θ::AbstractVector;
    opts=SolverOptions())
z_next = contact_step!(mechanism, z, θ; opts)
jacobian_state, jacobian_contact = Dojo.get_contact_gradients(mechanism)
return z_next, jacobian_state, jacobian_contact
end

function contact_step!(mechanism::Mechanism, z::AbstractVector, θ::AbstractVector;
    opts=SolverOptions())
set_data!(mechanism.contacts, θ)
nu = input_dimension(mechanism)
step!(mechanism, z, zeros(nu), opts=opts)
end

function loss(mechanism::Mechanism, θ::AbstractVector{T}, storage::Storage{T,N},
    timesteps::UnitRange{Int}; opts=SolverOptions(), derivatives::Bool=false) where {T,N}

nz = maximal_dimension(mechanism, attjac=true)
nd = length(θ)

Q = Diagonal([ones(3); 1e-1ones(3); ones(4); 1e-1ones(3)])
cost = 0.0
gradient = zeros(nd)
hessian = zeros(nd,nd)

d_contact = zeros(nz,nd)
z_prediction = get_maximal_state(storage, timesteps[1])
for i in timesteps
    z_true = get_maximal_state(storage, i+1)
    
    if derivatives
        z_prediction, ∂_state, ∂_contact = get_contact_gradients!(mechanism, z_prediction, θ; opts)
        attjac = attitude_jacobian(z_prediction, 1)
        d_contact = ∂_contact + ∂_state * d_contact
        gradient += (attjac * d_contact)' * Q * (z_prediction - z_true)
        hessian += (attjac * d_contact)' * Q * (attjac * d_contact)
    else
        z_prediction = contact_step!(mechanism, z_prediction, θ; opts)
    end
    cost += 0.5 * (z_prediction - z_true)'* Q *(z_prediction - z_true)
end
if derivatives
    return cost, gradient, hessian
else
    return cost
end
end