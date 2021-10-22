

################################################################################
using Test
using Random
using ForwardDiff
using Symbolics
using LinearAlgebra

function VRᵀmat(q::AbstractVector)
    SA[
        -q[2]  q[1] -q[4]  q[3];
        -q[3]  q[4]  q[1] -q[2];
        -q[4] -q[3]  q[2]  q[1];
    ]
end

function LVᵀmat(q::AbstractVector)
    SA[
        -q[2] -q[3] -q[4];
         q[1] -q[4]  q[3];
         q[4]  q[1] -q[2];
        -q[3]  q[2]  q[1];
    ]
end


function G(q::AbstractVector)
    qs = q[1]
    qv = q[2:4]
    m = zeros(eltype(q), 4, 3)
    m[1,:] += -qv
    m[2:end,:] += qs * I + skew(qv)
    return m
end

function ∂Rq_p(q::AbstractVector, p::AbstractVector)
    function f(q)
        return VRᵀmat(q) * LVᵀmat(q) * p
    end
    ∇ = ForwardDiff.jacobian(f, q)
    return ∇
end

Random.seed!(1)
p = [1, 2, 3.0]
q = [1., 3, 4, 2]
q ./= norm(q)
qq = UnitQuaternion(q...)

@variables qv[1:4], pv[1:3]
qu = [qv[1], qv[2], qv[3], qv[4]]
pu = [pv[1], pv[2], pv[3]]
y = VRᵀmat(qu) * LVᵀmat(qu) * pu
yjac = Symbolics.jacobian(y, qu)
yjac_func = eval(Symbolics.build_function(yjac, qu, pu)[1])
yjac[1,1]
yjac[1,2]

# Rotation
@test norm(qq * p - VRᵀmat(q) * LVᵀmat(q) * p) < 1e-12
# Jacobian
@test norm(2.0 * VRᵀmat(q) * LVᵀmat(q) * skew(-p) - ∂Rq_p(q, p) * G(q)) < 1e-12
@test norm(2.0 * VRᵀmat(q) * LVᵀmat(q) * skew(-p) - ∂vrotate∂q(p, UnitQuaternion(q...)) * G(q)) < 1e-12

# Attitude Jacobian
@test norm(LVᵀmat(q) - G(q)) < 1e-12

# yjac_func
@test norm(2.0 * VRᵀmat(q) * LVᵀmat(q) * skew(-p) - yjac_func(q, p) * LVᵀmat(q)) < 1e-12

# ∂Rq_p
@test norm(2.0 * VRᵀmat(q) * LVᵀmat(q) * skew(-p) - ∂Rq_p(q, p) * LVᵀmat(q)) < 1e-12

# # BQ
# @test norm(VRᵀmat(q) * LVᵀmat(q) * skew(-p) - BQ(q, p) * LVᵀmat(q)) < 1e-12
# @test norm(2 * BQ(q, p) - yjac_func(q, p)) < 1e-12

function B(z, p) 
    nx = 3 
    nq = nx
    q = z[4:7]
    Bxmat = SA[
        1 0 0
        0 1 0]
    [[szeros(1,nx); Bxmat] [szeros(1,nq); Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)*2.0]]
end

@variables z[1:7], p[1:3], b[1:3]
B(z, p)
Bb = transpose(B(z, p)) * b
dBb = Symbolics.jacobian(Bb, z)

println(dBb[4,4])
println(dBb[4,5])
println(dBb[4,6])
println(dBb[4,7])

println(dBb[5,4])
println(dBb[5,5])
println(dBb[5,6])
println(dBb[5,7])

println(dBb[6,4])
println(dBb[6,5])
println(dBb[6,6])
println(dBb[6,7])

function dB(x, q, b, p) 
    p₁, p₂, p₃ = p 
    z₁, z₂, z₃ = x
    z₄, z₅, z₆, z₇ = q 
    b₁, b₂, b₃ = b
    dBb = zeros(6, 7) 

    dBb[4,4] = b₂*(4.0p₂*z₆ + 4.0p₃*z₇) + b₃*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dBb[4,5] = b₂*(4.0p₂*z₇ - (4.0p₃*z₆)) + b₃*(4.0p₃*z₅ - (4.0p₂*z₄))
    dBb[4,6] = b₂*(4.0p₂*z₄ - (4.0p₃*z₅)) + b₃*(4.0p₂*z₇ - (4.0p₃*z₆))
    dBb[4,7] = b₂*(4.0p₂*z₅ + 4.0p₃*z₄) + b₃*(4.0p₂*z₆ + 4.0p₃*z₇)

    dBb[5,4] = b₂*(4.0p₃*z₄ - (4.0p₁*z₆)) + b₃*(4.0p₁*z₅ + 4.0p₃*z₇)
    dBb[5,5] = b₂*(4.0p₃*z₅ - (4.0p₁*z₇)) + b₃*(4.0p₁*z₄ + 4.0p₃*z₆)
    dBb[5,6] = b₂*(-4.0p₁*z₄ - (4.0p₃*z₆)) + b₃*(4.0p₃*z₅ - (4.0p₁*z₇))
    dBb[5,7] = b₂*(-4.0p₁*z₅ - (4.0p₃*z₇)) + b₃*(4.0p₃*z₄ - (4.0p₁*z₆))

    dBb[6,4] = b₂*(-4.0p₁*z₇ - (4.0p₂*z₄)) + b₃*(4.0p₁*z₄ - (4.0p₂*z₇))
    dBb[6,5] = b₂*(4.0p₁*z₆ - (4.0p₂*z₅)) + b₃*(-4.0p₁*z₅ - (4.0p₂*z₆))
    dBb[6,6] = b₂*(4.0p₁*z₅ + 4.0p₂*z₆) + b₃*(4.0p₁*z₆ - (4.0p₂*z₅))
    dBb[6,7] = b₂*(4.0p₂*z₇ - (4.0p₁*z₄)) + b₃*(-4.0p₁*z₇ - (4.0p₂*z₄))

    return dBb * [I zeros(3,3); zeros(4,3) G(q)]
end

x0 = rand(3) 
q0 = rand(4)
b0 = rand(3) 
p0 = rand(3) 
z0 = [x0; q0]
norm(dB(x0, q0, b0, p0) - ForwardDiff.jacobian(z -> transpose(B(z, p0)) * b0, z0) * [I zeros(3,3); zeros(4,3) G(q0)]) < 1.0e-8

function N(z, p) 
    q = z[4:7]
    ainv3 = [0.0; 0.0; 1.0]
    [transpose(ainv3) transpose(ainv3) * VRᵀmat(q)*LVᵀmat(q)*skew(-p)*2.0]
end

@variables γ[1:1]

Nγ = transpose(N(z, p)) * γ[1]

dNγ = Symbolics.jacobian(Nγ, z)

println(dNγ[4,4])
println(dNγ[4,5])
println(dNγ[4,6])
println(dNγ[4,7])

println(dNγ[5,4])
println(dNγ[5,5])
println(dNγ[5,6])
println(dNγ[5,7])

println(dNγ[6,4])
println(dNγ[6,5])
println(dNγ[6,6])
println(dNγ[6,7])

function dN(x, q, γ, p) 
    p₁, p₂, p₃ = p 
    z₁, z₂, z₃ = x
    γ₁ = γ[1]
    z₄, z₅, z₆, z₇ = q 
    dNγ = zeros(6, 7) 

    dNγ[4,4] = γ₁*(4.0p₂*z₄ - (4.0p₃*z₅))
    dNγ[4,5] = γ₁*(-4.0p₂*z₅ - (4.0p₃*z₄))
    dNγ[4,6] = γ₁*(-4.0p₂*z₆ - (4.0p₃*z₇))
    dNγ[4,7] = γ₁*(4.0p₂*z₇ - (4.0p₃*z₆))

    dNγ[5,4] = γ₁*(-4.0p₁*z₄ - (4.0p₃*z₆))
    dNγ[5,5] = γ₁*(4.0p₁*z₅ + 4.0p₃*z₇)
    dNγ[5,6] = γ₁*(4.0p₁*z₆ - (4.0p₃*z₄))
    dNγ[5,7] = γ₁*(4.0p₃*z₅ - (4.0p₁*z₇))

    dNγ[6,4] = γ₁*(4.0p₁*z₅ + 4.0p₂*z₆)
    dNγ[6,5] = γ₁*(4.0p₁*z₄ - (4.0p₂*z₇))
    dNγ[6,6] = γ₁*(4.0p₁*z₇ + 4.0p₂*z₄)
    dNγ[6,7] = γ₁*(4.0p₁*z₆ - (4.0p₂*z₅))

    return dNγ * [I zeros(3,3); zeros(4,3) G(q)]
end

x0 = rand(3) 
q0 = rand(4)
γ0 = rand(1) 
p0 = rand(3) 
z0 = [x0; q0]
norm(dN(x0, q0, γ0, p0) - ForwardDiff.jacobian(z -> transpose(N(z, p0)) * γ0[1], z0) * [I zeros(3,3); zeros(4,3) G(q0)]) < 1.0e-8

# 
@variables q[1:4] ω[1:3] p[1:3]
function Bω(q, ω, p) 
      Bxmat = SA[
                 1 0 0
                 0 1 0
                ]
    Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)*2.0 * ω
end

Bω_sym = Bω(q, ω, p)

dBω_sym = Symbolics.jacobian(Bω_sym, q) #* G(q)
dBω_sym = simplify.(dBω_sym)

println(dBω_sym[1, 1])
println(dBω_sym[2, 1])

println(dBω_sym[1, 2])
println(dBω_sym[2, 2])

println(dBω_sym[1, 3])
println(dBω_sym[2, 3])

println(dBω_sym[1, 4])
println(dBω_sym[2, 4])

function dBω(q, ω, p)
    q₁, q₂, q₃, q₄ = q 
    ω₁, ω₂, ω₃ = ω 
    p₁, p₂, p₃ = p 
    dB = zeros(2, 4)

    dB[1, 1] = ω₁*(4.0p₂*q₃ + 4.0p₃*q₄) + ω₂*(4.0p₃*q₁ - (4.0p₁*q₃)) + ω₃*(-4.0p₁*q₄ - (4.0p₂*q₁))
    dB[2, 1] = ω₁*(-4.0p₂*q₂ - (4.0p₃*q₁)) + ω₂*(4.0p₁*q₂ + 4.0p₃*q₄) + ω₃*(4.0p₁*q₁ - (4.0p₂*q₄))

    dB[1, 2] = ω₁*(4.0p₂*q₄ - (4.0p₃*q₃)) + ω₂*(4.0p₃*q₂ - (4.0p₁*q₄)) + ω₃*(4.0p₁*q₃ - (4.0p₂*q₂))
    dB[2, 2] = ω₁*(4.0p₃*q₂ - (4.0p₂*q₁)) + ω₂*(4.0p₁*q₁ + 4.0p₃*q₃) + ω₃*(-4.0p₁*q₂ - (4.0p₂*q₃))

    dB[1, 3] = ω₁*(4.0p₂*q₁ - (4.0p₃*q₂)) + ω₂*(-4.0p₁*q₁ - (4.0p₃*q₃)) + ω₃*(4.0p₁*q₂ + 4.0p₂*q₃)
    dB[2, 3] = ω₁*(4.0p₂*q₄ - (4.0p₃*q₃)) + ω₂*(4.0p₃*q₂ - (4.0p₁*q₄)) + ω₃*(4.0p₁*q₃ - (4.0p₂*q₂))

    dB[1, 4] = ω₁*(4.0p₂*q₂ + 4.0p₃*q₁) + ω₂*(-4.0p₁*q₂ - (4.0p₃*q₄)) + ω₃*(4.0p₂*q₄ - (4.0p₁*q₁))
    dB[2, 4] = ω₁*(4.0p₂*q₃ + 4.0p₃*q₄) + ω₂*(4.0p₃*q₁ - (4.0p₁*q₃)) + ω₃*(-4.0p₁*q₄ - (4.0p₂*q₁))

    return dB
end

q0 = rand(4) 
q0 ./= norm(q0)
ω0 = rand(3) 
p0 = rand(3) 

@test norm(dBω(q0, ω0, p0) - ForwardDiff.jacobian(a -> Bω(a, ω0, p0), q0)) < 1.0e-8
