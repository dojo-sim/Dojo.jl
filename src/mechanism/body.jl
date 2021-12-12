mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String

    m::T
    J::SMatrix{3,3,T,9}

    state::State{T}

    shape::Shape{T}

    parentid::Int64
    childid::Int64

    function Body(m::Real, J::AbstractArray; name::String="", shape::Shape=EmptyShape())
        T = promote_type(eltype.((m, J))...)
        new{T}(getGlobalID(), name, m, J, State{T}(), shape, 0, 0)
    end

    function Body{T}(contents...) where T
        new{T}(getGlobalID(), contents...)
    end
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64
    name::String

    shape::Shape{T}


    function Origin{T}(; name::String="", shape::Shape=EmptyShape()) where T
        new{T}(getGlobalID(), name, shape)
    end

    function Origin(; name::String="", shape::Shape=EmptyShape())
        Origin{Float64}(; name = name, shape = shape)
    end

    function Origin(body::Body{T}) where T
        new{T}(body.id, body.name, body.shape)
    end
end

function Base.deepcopy(b::Body{T}) where T
    contents = []
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    return Body{T}(contents...)
end

Base.length(::Body) = 6
Base.zero(::Body{T}) where T = szeros(T, 6, 6)

@inline getM(body::Body{T}) where T = [[I * body.m; szeros(T,3,3)] [szeros(T,3,3); body.J]]

# Derivatives for linearizations
function ∂F∂z(body::Body{T}, Δt::T; attjac::Bool = true) where T
    state = body.state
    q2 = state.q2[1]
    ϕ25 = state.ϕsol[2]
    Z3 = szeros(T,3,3)
    Z34 = szeros(T,3,4)
    ZT = attjac ? szeros(T,6,6) : szeros(T,6,7)
    ZR = szeros(T,7,6)

    x1, q1 = posargs1(state)
    x2, q2 = posargs2(state)
    x3, q3 = posargs3(state, Δt)

    AposT = [-I Z3]
    AvelT = [Z3 -I*body.m] # solving for impulses

    AposR = [-∂integrator∂q(q2, ϕ25, Δt, attjac = attjac) szeros(4,3)]

    J = body.J
    
    rot_q1(q) = 2/Δt * LVᵀmat(UnitQuaternion(q..., false))' * Tmat() * Rmat(q2)' * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(q2)
    rot_q2(q) = 2/Δt * LVᵀmat(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt))' * Lmat(UnitQuaternion(q..., false)) * Vᵀmat() * body.J * Vmat() * Lmat(UnitQuaternion(q..., false))' * vector(getq3(UnitQuaternion(q..., false), state.ϕsol[2], Δt)) + 2/Δt * LVᵀmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt))' * Tmat() * Rmat(UnitQuaternion(q..., false))' * Vᵀmat() * body.J * Vmat() * Lmat(getq3(UnitQuaternion(q..., false), -state.ϕ15, Δt))' * q
    dynR_ϕ15 = -1.0 * FiniteDiff.finite_difference_jacobian(rot_q1, vector(q1)) * ∂integrator∂ϕ(q2, -state.ϕ15, Δt)
    dynR_q2 = FiniteDiff.finite_difference_jacobian(rot_q2, vector(q2))
    AvelR = attjac ? [dynR_q2 * LVᵀmat(q2) dynR_ϕ15] : [dynR_q2 dynR_ϕ15]
    
    return [[AposT;AvelT] ZT;
             ZR [AposR;AvelR]]
end

function ∂F∂u(body::Body{T}, Δt) where T
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)

    BposT = [Z3 Z3]
    BvelT = [-I Z3]
    BposR = [Z43 Z43]
    BvelR = [Z3 -I]
    return [BposT;BvelT;BposR;BvelR]
end
