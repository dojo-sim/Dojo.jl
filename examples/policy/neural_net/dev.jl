using Test
using LinearAlgebra
using Symbolics
using BenchmarkTools

include("layer.jl")
include("net.jl")


mutable struct Net147{N,NI,Nθ,NO,T}
    dim_layers::Vector{Int}
    dim_parameters::Vector{Int}
    layers::Vector{<:Layer147{T}}
    θ::Vector{T}
    ∂xo∂xi_stages::Vector{Matrix{T}}
    ∂xo∂θ_stages::Vector{Matrix{T}}
    ∂xo∂xi::Matrix{T}
    ∂xo∂θ::Matrix{T}
end

function Net147(NI, NO; dim_layers=zeros(Int,0), activations=[x->x for i=1:length(dim_layers)+1], T=Float64)
    @assert length(dim_layers) == length(activations) - 1
    N = 1 + length(dim_layers)
    dim_layers = [dim_layers; NO]

    layers = Vector()
    ∂xo∂xi_stages = Vector{Matrix{T}}()
    ∂xo∂θ_stages = Vector{Matrix{T}}()
    ni = NI
    for i = 1:N
        no = dim_layers[i]
        layer = Layer147(ni, no, activations[i])
        push!(layers, layer)
        push!(∂xo∂xi_stages, zeros(NO, ni))
        push!(∂xo∂θ_stages, zeros(NO, layer.nθ))
        ni = no
    end
    layers = [layers...]
    dim_parameters = [layer.nθ for layer in layers]
    Nθ = sum(dim_parameters)

    θ = zeros(Nθ)
    ∂xo∂xi = zeros(NO, NI)
    ∂xo∂θ = zeros(NO, Nθ)

    T = eltype(layers[1].xi)
    return Net147{N,NI,Nθ,NO,T}(dim_layers, dim_parameters, layers,
        θ,
        ∂xo∂xi_stages, ∂xo∂θ_stages, ∂xo∂xi, ∂xo∂θ)
end

parameter_dimension(net::Net147{N,NI,Nθ}) where {N,NI,Nθ} = Nθ
get_input(net::Net147) = net.layers[1].xi
get_output(net::Net147{N}) where N = net.layers[N].xo

function set_input!(net::Net147, xi::Vector{T}) where T
    net.layers[1].xi = xi
    return nothing
end

function set_parameters!(net::Net147{N}, θ::Vector{T}) where {N,T}
    net.θ = θ
    layers = net.layers
    off = 0
    for i = 1:N
        nθ = layers[i].nθ
        for j = 1:nθ
            layers[i].θ[j] = θ[off + j]
        end
        off += nθ
    end
    return nothing
end

function evaluation!(net::Net147, xi::Vector{T}, θ::Vector{T}) where T
    set_input!(net, xi)
    set_parameters!(net, θ)
    evaluation!(net)
    return nothing
end

function evaluation!(net::Net147{N}) where N
    for i = 1:N
        evaluation!(net.layers[i])
        (i < N) && (net.layers[i+1].xi = net.layers[i].xo)
    end
    return nothing
end

function jacobian_input!(net::Net147{N}) where N
    # evaluation
    for i = 1:N
        layer = net.layers[i]
        evaluation!(layer)
        jacobian_input!(layer)
        (i < N) && (net.layers[i+1].xi = layer.xo)
    end

    # jacobian
    layer = net.layers[N]
    net.∂xo∂xi_stages[N] = layer.∂xo∂xi
    for i = N-1:-1:1
        layer = net.layers[i]
        mul!(net.∂xo∂xi_stages[i], net.∂xo∂xi_stages[i+1], layer.∂xo∂xi)
    end
    net.∂xo∂xi = net.∂xo∂xi_stages[1]
    return nothing
end

function jacobian_input!(net::Net147, xi::Vector{T}, θ::Vector{T}) where T
    set_input!(net, xi)
    set_parameters!(net, θ)
    jacobian_input!(net)
    return nothing
end

function jacobian_parameters!(net::Net147{N}) where N
    # evaluation & jacobian input
    jacobian_input!(net)

    # jacobian parameters
    for i = 1:N
        layer = net.layers[i]
        jacobian_parameters!(layer)
    end
    for i = 1:N-1
        layer = net.layers[i]
        mul!(net.∂xo∂θ_stages[i], net.∂xo∂xi_stages[i+1], layer.∂xo∂θ)
    end
    net.∂xo∂θ_stages[N] = net.layers[N].∂xo∂θ

    off = 0
    for i = 1:N
        nθ = net.layers[i].nθ
        net.∂xo∂θ[:, off .+ (1:nθ)] = net.∂xo∂θ_stages[i]
        off += nθ
    end
    return nothing
end

function jacobian_parameters!(net::Net147, xi::Vector{T}, θ::Vector{T}) where T
    set_input!(net, xi)
    set_parameters!(net, θ)
    jacobian_parameters!(net)
    return nothing
end


ni = 16
no = 18
xi0 = rand(ni)
net0 = Net147(ni, no, dim_layers=[5,5], activations=[x->tanh.(x), x->tanh.(x), x->x])
θ0 = rand(parameter_dimension(net0))

jacobian_input!(net0, xi0, θ0)
jacobian_parameters!(net0, xi0, θ0)
@benchmark $jacobian_parameters!($net0, $xi0, $θ0)


jacobian_parameters!(net0, xi0, θ0)
net0.∂xo∂θ

function local_evaluation(net, xi, θ)
    evaluation!(net, xi, θ)
    return copy(get_output(net))
end

J0 = FiniteDiff.finite_difference_jacobian(θ0 -> local_evaluation(net0, xi0, θ0), θ0)
@test norm(net0.∂xo∂θ - J0, Inf) < 1e-3



net0.∂xo∂θ
J0
plot(Gray.(abs.(net0.∂xo∂θ)))
plot(Gray.(abs.(J0)))


jacobian_input!(net0, xi0, θ0)
net0.∂xo∂xi

function local_evaluation(net, xi, θ)
    evaluation!(net, xi, θ)
    return copy(get_output(net))
end

J0 = FiniteDiff.finite_difference_jacobian(xi0 -> local_evaluation(net0, xi0, θ0), xi0)
@test norm(net0.∂xo∂xi - J0) < 1e-3



# @benchmark $jacobian_input!($net0)





# evaluation!(net0, 12xi0, θ0)
evaluation!(net0, 0.1*ones(3), ones(285))
net0.layers[1].xi
net0.layers[1].θ
net0.layers[2].θ
net0.layers[3].θ
net0.layers[4].θ
net0.layers[1].xo
get_output(net0)

@benchmark $evaluation!($net0, $xi0, $θ0)




set_parameters!(net0, θ0)
Main.@code_warntype set_parameters!(net0, θ0)
@benchmark $set_parameters!($net0, $θ0)

set_input!(net0, xi0)
@benchmark $set_input!($net0, $xi0)
evaluation!(net0)
@benchmark $evaluation!($net0)


N = 100
inputs = [rand(3) for i = 1:N]
outputs = [[input; sin.(input)] for input in inputs]

layer1 = Layer147(3, 6, x -> x)
slow_evaluation(inputs[1], layer1)
evaluation!(layer1)
jacobian_input!(layer1)
Main.@code_warntype jacobian_input!(layer1)

@benchmark $evaluation!($layer1)
@benchmark $jacobian_input!($layer1)
@benchmark $jacobian_parameters!($layer1)
