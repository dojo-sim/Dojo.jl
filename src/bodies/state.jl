"""
    State{T} 

    state information in maximal coordinates for Body at time steps: 1, 2, 3. 
    information at time step 3 is recovered using configurations at time step 2 and velocities at time step 2.5.

    x1: position at previous time step
    q1: orientation (Quaternion) at previous time step 
    v15: linear velocity at time step 1.5 (midpoint)
    ϕ15: angular velocity at time step 1.5 (midpoint)

    x2: position at current time step 
    q2: orientation (Quaternion) at current time step 
    F2: linear impulse (force * time step) applied at current time step 
    τ2: angular impulse (torque * timestep) applied at current time step 

    vsol: linear velocity at time step 2.5 (midpoint); contains current value (index 1) and candidate value (index 2)
    ϕsol: angular velocity at time step 2.5 (midpoint); contains current value (index 1) and candidate value (index 2)

    d: implicit dynamics evaluator (zero vector indicates physics are satisfied)
    D: Jacobian of implicit dynamics
"""
mutable struct State{T}
    # previous state
    x1::SVector{3,T}
    q1::Quaternion{T}
    v15::SVector{3,T}
    ϕ15::SVector{3,T}

    # current state
    x2::SVector{3,T}
    q2::Quaternion{T}
    F2::SVector{3,T}
    τ2::SVector{3,T}

    # solution estimate [before step; after step]
    vsol::Vector{SVector{3,T}}
    ϕsol::Vector{SVector{3,T}}

    # dynamics and Jacobian evaluation
    d::SVector{6,T}
    D::SMatrix{6,6,T,36}

    function State{T}() where T
        x1 = szeros(T, 3)
        q1 = one(Quaternion{T})
        v15 = szeros(T, 3)
        ϕ15 = szeros(T, 3)

        x2 = szeros(T, 3)
        q2 = one(Quaternion{T})
        F2 = szeros(T, 3)
        τ2 = szeros(T, 3)

        vsol = [szeros(T, 3) for i=1:2]
        ϕsol = [szeros(T, 3) for i=1:2]

        d = szeros(T, 6)
        D = szeros(T, 6, 6)

        new{T}(x1, q1, v15, ϕ15, x2, q2, F2, τ2, vsol, ϕsol, d, D)
    end
end
