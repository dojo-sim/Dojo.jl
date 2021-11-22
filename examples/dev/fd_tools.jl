#
# function fd_diagonal∂damper∂ʳvel(j::Joint{T}, x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector,
#     x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1a)
#         v1a = vω1a[1:3]
#         ω1a = vω1a[4:6]
#         Fτa = damperforcea(j, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b)
#         return Fτa
#     end
#     return fdjac(vω1a -> f(vω1a), Vector([v1a; ω1a]))
# end
#
# function fd_offdiagonal∂damper∂ʳvel(j::Joint{T}, x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector,
#     x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1b)
#         v1b = vω1b[1:3]
#         ω1b = vω1b[4:6]
#         Fτb = damperforcea(j, x1a, v1a, q1a, ω1a, x1b, v1b, q1b, ω1b)
#         return Fτb
#     end
#     return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
# end
#
# function fd_offdiagonal∂damper∂ʳvel(j::Joint{T}, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1b)
#         v1b = vω1b[1:3]
#         ω1b = vω1b[4:6]
#         Fτb = -1.0 * damperforceb(j, x1b, v1b, q1b, ω1b)
#         return Fτb
#     end
#     return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
# end
#
#
# function fd_diagonal∂spring∂ʳvel(j::Joint{T}, x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector,
#     x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1a)
#         v1a = vω1a[1:3]
#         ω1a = vω1a[4:6]
#         Fτa = springforcea(j, x1a, q1a, x1b, q1b)
#         return Fτa
#     end
#     return fdjac(vω1a -> f(vω1a), Vector([v1a; ω1a]))
# end
#
# function fd_offdiagonal∂spring∂ʳvel(j::Joint{T}, x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector,
#     x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1b)
#         v1b = vω1b[1:3]
#         ω1b = vω1b[4:6]
#         Fτb = springforceb(j, x1a, q1a, x1b, q1b)
#         return Fτb
#     end
#     return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
# end
#
# function fd_offdiagonal∂spring∂ʳvel(j::Joint{T}, x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector) where T
#     function f(vω1b)
#         v1b = vω1b[1:3]
#         ω1b = vω1b[4:6]
#         Fτb = springforceb(j, x1b, q1b)
#         return Fτb
#     end
#     return fdjac(vω1b -> f(vω1b), Vector([v1b; ω1b]))
# end
#
#
#
# ## Spring velocity derivatives
# @inline function diagonal∂spring∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     # A = nullspacemat(joint)
#     # AᵀA = zerodimstaticadjoint(A) * A
#     # X, Q = ∂g∂ʳposa(joint, xa, qa, xb, qb)
#     # Fv = AᵀA * joint.spring * AᵀA * X
#     # Fω = AᵀA * joint.spring * AᵀA * Q
#     Z = szeros(T, 3, 3)
#     # return [[Fv; Z] [Fω; Z]]
#     return [[Z; Z] [Z; Z]]
# end
# @inline function offdiagonal∂spring∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     # A = nullspacemat(joint)
#     # AᵀA = zerodimstaticadjoint(A) * A
#     Z = szeros(T, 3, 3)
#     # return [[AᵀA * joint.spring * AᵀA; Z] [Z; Z]]
#     return [[Z; Z] [Z; Z]]
# end
# @inline function offdiagonal∂spring∂ʳvel(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
#     # A = nullspacemat(joint)
#     # AᵀA = zerodimstaticadjoint(A) * A
#     Z = szeros(T, 3, 3)
#     # return [[AᵀA * joint.spring * AᵀA; Z] [Z; Z]]
#     return [[Z; Z] [Z; Z]]
# end
#
# ## Damper velocity derivatives
# @inline function diagonal∂damper∂ʳvel(joint::Rotational{T}) where T # never used
#     A = nullspacemat(joint)
#     AᵀA = zerodimstaticadjoint(A) * A
#     Z = szeros(T, 3, 3)
#     return [[Z; Z] [Z; -2 * AᵀA * joint.damper * AᵀA]]
# end
# @inline function offdiagonal∂damper∂ʳvel(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     # invqbqa = q2b\q2a
#     # A = nullspacemat(joint)
#     # AᵀA = zerodimstaticadjoint(A) * A
#     # Z = szeros(T, 3, 3)
#     # return [[Z; Z] [Z; 2*VLmat(invqbqa)*RVᵀmat(invqbqa)* AᵀA * Diagonal(joint.damper) * AᵀA]]
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     C = 2 * Aᵀ * A * joint.damper * Aᵀ * A
#     δq = qa \ qb
#     # invδq = qb \ qa
#     Z = szeros(T, 3, 3)
#     # return [[Z; Z] [Z; VRᵀmat(invδq) * LVᵀmat(invδq) * C * VRᵀmat(δq) * LVᵀmat(δq)]]
#     return [[Z; Z] [Z; C * VRᵀmat(δq) * LVᵀmat(δq)]]
#     # offdiagonal∂damper∂ʳvel(joint, xb, qa \ qb)
# end
# @inline function offdiagonal∂damper∂ʳvel(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
#     # invqb = inv(q2b)
#     # A = nullspacemat(joint)
#     # AᵀA = zerodimstaticadjoint(A) * A
#     # Z = szeros(T, 3, 3)
#     # return [[Z; Z] [Z; 2*VLmat(invqb)*RVᵀmat(invqb)* AᵀA * Diagonal(joint.damper) * AᵀA]]
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     C = 2 * Aᵀ * A * joint.damper * Aᵀ * A
#     q = qb
#     invq = inv(qb)
#     Z = szeros(T, 3, 3)
#     return [[Z; Z] [Z; VRᵀmat(invq) * LVᵀmat(invq) * C * VRᵀmat(q) * LVᵀmat(q)]]
#     # offdiagonal∂damper∂ʳvel(joint, zeros(3), UnitQuaternion(1,0,0,0.0), xb, qb)
# end
#
#
#
#
#
# ## Damper velocity derivatives
# # @inline function data_diagonal∂damper∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
# #     A = nullspacemat(joint)
# #     AᵀA = zerodimstaticadjoint(A) * A
# #     C = 2 * AᵀA * joint.damper * AᵀA
# #     Fq = C * ∂vrot∂q(ωb, qa \ qb) *
#
# #     Z = szeros(T, 3, 3)
# #     return [[Z; Z] [Fq; Z]]
# # end
# @inline function data_offdiagonal∂damper∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
#     Z = szeros(T, 3, 3)
#     return [[Z; Z] [Z; Z]]
# end
# @inline function data_offdiagonal∂damper∂ʳvel(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
#     Z = szeros(T, 3, 3)
#     return [[Z; Z] [Z; Z]]
# end
#
#
#
#
#
#
#
#
# @inline function damperforcea(joint::Rotational, xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     velocity = A * (vrotate(ωb,qa\qb) - ωa) # in body1's frame
#     force = 2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
#     return [szeros(3);force]
# end
# @inline function damperforceb(joint::Rotational, xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#     xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#
#     velocity = A * (vrotate(ωb,qa\qb) - ωa) # in body1's frame
#
#     force = -2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
#     force = vrotate(force,qb\qa) # in body2's frame
#     return [szeros(3);force]
# end
# @inline function damperforceb(joint::Rotational, xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#
#     velocity = A * vrotate(ωb,qb)  # in world frame
#
#     force = -2 * Aᵀ * A * joint.damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
#     force = vrotate(force,inv(qb)) # in body2's frame
#     return [szeros(3);force]
# end
#
#
# # ## Damper velocity derivatives
# # @inline function data_diagonal∂damper∂ʳvel(joint::Rotational{T}) where T # never used
# #     A = nullspacemat(joint)
# #     AᵀA = zerodimstaticadjoint(A) * A
# #     Z = szeros(T, 3, 3)
# #     Fqa = C * ...
# #     C = 2 * AᵀA joint.damper * AᵀA *
# #     return [[Z; Z] [Z; Fqa]]
# # end
#
# # @inline function data_offdiagonal∂damper∂ʳvel(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
# #     # invqbqa = q2b\q2a
# #     # A = nullspacemat(joint)
# #     # AᵀA = zerodimstaticadjoint(A) * A
# #     # Z = szeros(T, 3, 3)
# #     # return [[Z; Z] [Z; 2*VLmat(invqbqa)*RVᵀmat(invqbqa)* AᵀA * Diagonal(joint.damper) * AᵀA]]
# #     A = nullspacemat(joint)
# #     Aᵀ = zerodimstaticadjoint(A)
# #     C = 2 * Aᵀ * A * joint.damper * Aᵀ * A
# #     δq = qa \ qb
# #     # invδq = qb \ qa
# #     Z = szeros(T, 3, 3)
# #     # return [[Z; Z] [Z; VRᵀmat(invδq) * LVᵀmat(invδq) * C * VRᵀmat(δq) * LVᵀmat(δq)]]
# #     return [[Z; Z] [Z; C * VRᵀmat(δq) * LVᵀmat(δq)]]
# #     # offdiagonal∂damper∂ʳvel(joint, xb, qa \ qb)
# # end
# # @inline function data_offdiagonal∂damper∂ʳvel(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
# #     # invqb = inv(q2b)
# #     # A = nullspacemat(joint)
# #     # AᵀA = zerodimstaticadjoint(A) * A
# #     # Z = szeros(T, 3, 3)
# #     # return [[Z; Z] [Z; 2*VLmat(invqb)*RVᵀmat(invqb)* AᵀA * Diagonal(joint.damper) * AᵀA]]
# #     A = nullspacemat(joint)
# #     Aᵀ = zerodimstaticadjoint(A)
# #     C = 2 * Aᵀ * A * joint.damper * Aᵀ * A
# #     q = qb
# #     invq = inv(qb)
# #     Z = szeros(T, 3, 3)
# #     return [[Z; Z] [Z; VRᵀmat(invq) * LVᵀmat(invq) * C * VRᵀmat(q) * LVᵀmat(q)]]
# #     # offdiagonal∂damper∂ʳvel(joint, zeros(3), UnitQuaternion(1,0,0,0.0), xb, qb)
# # end
#
#
#
# #
# # ## Spring velocity derivatives
# # @inline function data_diagonal∂spring∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
# #     A = nullspacemat(joint)
# #     AᵀA = zerodimstaticadjoint(A) * A
# #     X, Q = ∂g∂ʳposa(joint, xa, qa, xb, qb)
# #     Fv = AᵀA * joint.spring * AᵀA * X
# #     Fω = AᵀA * joint.spring * AᵀA * Q
# #     Z = szeros(T, 3, 3)
# #     return [[Fv; Z] [Fω; Z]]
# #     # return [[Z; Z] [Z; Z]]
# # end
# # @inline function data_offdiagonal∂spring∂ʳvel(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
# #     A = nullspacemat(joint)
# #     AᵀA = zerodimstaticadjoint(A) * A
# #     Z = szeros(T, 3, 3)
# #     return [[AᵀA * joint.spring * AᵀA; Z] [Z; Z]]
# #     # return [[Z; Z] [Z; Z]]
# # end
# # @inline function data_offdiagonal∂spring∂ʳvel(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
# #     A = nullspacemat(joint)
# #     AᵀA = zerodimstaticadjoint(A) * A
# #     Z = szeros(T, 3, 3)
# #     return [[AᵀA * joint.spring * AᵀA; Z] [Z; Z]]
# #     # return [[Z; Z] [Z; Z]]
# # end
