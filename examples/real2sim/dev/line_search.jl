# TODO: Implement safeguards


"""
`StrongWolfe`: This linesearch algorithm guarantees that the step length
satisfies the (strong) Wolfe conditions.
See Nocedal and Wright - Algorithms 3.5 and 3.6

This algorithm is mostly of theoretical interest, users should most likely
use `MoreThuente`, `HagerZhang` or `BackTracking`.

## Parameters:  (and defaults)
* `c_1 = 1e-4`: Armijo condition
* `c_2 = 0.9` : second (strong) Wolfe condition
* `ρ = 2.0` : bracket growth
"""
@with_kw struct StrongWolfe{T}
    c_1::T = 1e-4
    c_2::T = 0.9
    ρ::T = 2.0
end


function StrongWolfe(f, g, x::AbstractArray{T}, p::AbstractArray{T}, α::Real,
    x_new::AbstractArray{T}, ϕ_0, dϕ_0) where T

    ϕ, dϕ, ϕdϕ = make_ϕ_dϕ_ϕdϕ(f, g, x_new, x, p)
    ls(ϕ, dϕ, ϕdϕ, α, ϕ_0, dϕ_0)
end


function (ls::StrongWolfe)(ϕ, dϕ, ϕdϕ,
                           alpha0::T, ϕ_0, dϕ_0) where T
    @unpack c_1, c_2, ρ = ls

    zeroT = convert(T, 0)

    # Step-sizes
    a_0 = zeroT
    a_iminus1 = a_0
    a_i = alpha0
    a_max = convert(T, 65536)

    # ϕ(alpha) = df.f(x + alpha * p)
    ϕ_a_iminus1 = ϕ_0
    ϕ_a_i = convert(T, NaN)

    # ϕ'(alpha) = dot(g(x + alpha * p), p)
    dϕ_a_i = convert(T, NaN)

    # Iteration counter
    i = 1

    while a_i < a_max
        ϕ_a_i = ϕ(a_i)

        # Test Wolfe conditions
        if (ϕ_a_i > ϕ_0 + c_1 * a_i * dϕ_0) ||
            (ϕ_a_i >= ϕ_a_iminus1 && i > 1)
            a_star = zoom(a_iminus1, a_i,
                          dϕ_0, ϕ_0,
                          ϕ, dϕ, ϕdϕ)
            return a_star, ϕ(a_star)
        end

        dϕ_a_i = dϕ(a_i)

        # Check condition 2
        if abs(dϕ_a_i) <= -c_2 * dϕ_0
            return a_i, ϕ_a_i
        end

        # Check condition 3
        if dϕ_a_i >= zeroT # FIXME untested!
            a_star = zoom(a_i, a_iminus1,
                          dϕ_0, ϕ_0, ϕ, dϕ, ϕdϕ)
            return a_star, ϕ(a_star)
        end

        # Choose a_iplus1 from the interval (a_i, a_max)
        a_iminus1 = a_i
        a_i *= ρ

        # Update ϕ_a_iminus1
        ϕ_a_iminus1 = ϕ_a_i

        # Update iteration count
        i += 1
    end

    # Quasi-error response TODO make this error instead
    return a_max, ϕ(a_max)
end

function zoom(a_lo::T,
              a_hi::T,
              dϕ_0::Real,
              ϕ_0::Real,
              ϕ,
              dϕ,
              ϕdϕ,
              c_1::Real = convert(T, 1)/10^4,
              c_2::Real = convert(T, 9)/10) where T

    zeroT = convert(T, 0)
    # Step-size
    a_j = convert(T, NaN)

    # Count iterations
    iteration = 0
    max_iterations = 10

    # Shrink bracket
    while iteration < max_iterations
        iteration += 1

        ϕ_a_lo, ϕprime_a_lo = ϕdϕ(a_lo)

        ϕ_a_hi, ϕprime_a_hi = ϕdϕ(a_hi)

        # Interpolate a_j
        if a_lo < a_hi
            a_j = interpolate(a_lo, a_hi,
                              ϕ_a_lo, ϕ_a_hi,
                              ϕprime_a_lo, ϕprime_a_hi)
        else
            # TODO: Check if this is needed
            a_j = interpolate(a_hi, a_lo,
                              ϕ_a_hi, ϕ_a_lo,
                              ϕprime_a_hi, ϕprime_a_lo)
        end

        # Evaluate ϕ(a_j)
        ϕ_a_j = ϕ(a_j)

        # Check Armijo
        if (ϕ_a_j > ϕ_0 + c_1 * a_j * dϕ_0) ||
            (ϕ_a_j > ϕ_a_lo)
            a_hi = a_j
        else
            # Evaluate ϕprime(a_j)
            ϕprime_a_j = dϕ(a_j)

            if abs(ϕprime_a_j) <= -c_2 * dϕ_0
                return a_j
            end

            if ϕprime_a_j * (a_hi - a_lo) >= zeroT
                a_hi = a_lo
            end

            a_lo = a_j
        end
    end

    # Quasi-error response
    return a_j
end

# a_lo = a_{i - 1}
# a_hi = a_{i}
function interpolate(a_i1::Real, a_i::Real,
                     ϕ_a_i1::Real, ϕ_a_i::Real,
                     dϕ_a_i1::Real, dϕ_a_i::Real)
    d1 = dϕ_a_i1 + dϕ_a_i -
        3 * (ϕ_a_i1 - ϕ_a_i) / (a_i1 - a_i)
    d2 = sqrt(d1 * d1 - dϕ_a_i1 * dϕ_a_i)
    return a_i - (a_i - a_i1) *
        ((dϕ_a_i + d2 - d1) /
         (dϕ_a_i - dϕ_a_i1 + 2 * d2))
end


function make_ϕ_dϕ_ϕdϕ(f, g, x_new, x, s)
    function ϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Calculate ϕ'(a_i)
        f(x_new)
    end
    function dϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Calculate ϕ'(a_i)
        return g(x_new)' * s
    end
    function ϕdϕ(α)
        # Move a distance of alpha in the direction of s
        x_new .= x .+ α.*s

        # Calculate ϕ'(a_i)
        f(x_new), g(x_new)' * s
    end

    ϕ, dϕ, ϕdϕ
end
