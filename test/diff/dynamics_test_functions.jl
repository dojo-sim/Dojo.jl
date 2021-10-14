using ConstrainedDynamics
using ConstrainedDynamics: skewplusdiag



function dynTvel(vars)
    ezg = [0; 0; -9.81]
    Δt = 0.01
    m = 1.

    vc1 = vars[1:3]
    vc2 = vars[4:6]

    return m * ((vc2 - vc1) / Δt + ezg)
end

function dynRvel(vars)
    Δt = 0.01

    J = diagm(ones(3))
    ωc1 = vars[1:3]
    ωc2 = vars[4:6]
    sq2 = sqrt(4 / Δt^2 - ωc2' * ωc2)
    sq1 = sqrt(4 / Δt^2 - ωc1' * ωc1)
    
    return skewplusdiag(ωc2, sq2) * (J * ωc2) - skewplusdiag(ωc1, sq1) * (J * ωc1)
end