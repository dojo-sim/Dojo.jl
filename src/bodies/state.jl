"""
    State{T} 

    state information in maximal coordinates for Body at time steps: 1, 2, 3. 
    information at time step 3 is recovered using configurations at time step 2 and velocities at time step 2.5.

    x1: position at previous time step
    q1: orientation (Quaternion) at previous time step 
    v15: linear velocity at time step 1.5 (midpoint)
    ω15: angular velocity at time step 1.5 (midpoint)

    x2: position at current time step 
    q2: orientation (Quaternion) at current time step 
    JF2: linear impulse (force * time step) applied at current time step 
    Jτ2: angular impulse (torque * timestep) applied at current time step 

    vsol: linear velocity at time step 2.5 (midpoint); contains current value (index 1) and candidate value (index 2)
    ωsol: angular velocity at time step 2.5 (midpoint); contains current value (index 1) and candidate value (index 2)

    d: implicit dynamics evaluator (zero vector indicates physics are satisfied)
    D: Jacobian of implicit dynamics
"""
mutable struct State{T}
    # previous state
    x1::SVector{3,T}
    q1::Quaternion{T}
    v15::SVector{3,T}
    ω15::SVector{3,T}

    # current state
    x2::SVector{3,T}
    q2::Quaternion{T}
    JF2::SVector{3,T}
    Jτ2::SVector{3,T}

    # solution estimate [before step; after step]
    vsol::Vector{SVector{3,T}}
    ωsol::Vector{SVector{3,T}}

    # dynamics and Jacobian evaluation
    d::SVector{6,T}
    D::SMatrix{6,6,T,36}

    function State{T}() where T
        x1 = szeros(T, 3)
        q1 = one(Quaternion{T})
        v15 = szeros(T, 3)
        ω15 = szeros(T, 3)

        x2 = szeros(T, 3)
        q2 = one(Quaternion{T})
        JF2 = szeros(T, 3)
        Jτ2 = szeros(T, 3)

        vsol = [szeros(T, 3) for i=1:2]
        ωsol = [szeros(T, 3) for i=1:2]

        d = szeros(T, 6)
        D = szeros(T, 6, 6)

        new{T}(x1, q1, v15, ω15, x2, q2, JF2, Jτ2, vsol, ωsol, d, D)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, state::State{T}) where {T}
    summary(io, state)
    println(io, "")
    println(io, "x1:   "*string(state.x1))
    println(io, "q1:   "*string(state.q1))
    println(io, "v15:  "*string(state.v15))
    println(io, "ω15:  "*string(state.ω15))
    println(io, "x2:   "*string(state.x2))
    println(io, "q2:   "*string(state.q2))
    println(io, "JF2:   "*string(state.JF2))
    println(io, "Jτ2:   "*string(state.Jτ2))
    println(io, "vsol: "*string(state.vsol))
    println(io, "ωsol: "*string(state.ωsol))
    println(io, "d:    "*string(state.d))
    println(io, "D:    "*string(state.D))
end
