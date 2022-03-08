using MeshCat
using Colors

using EnhancedGJK 
using CoordinateTransformations
using StaticArrays

using GeometryTypes
using GeometryBasics

c1 = GeometryTypes.HyperRectangle(GeometryTypes.Vec(0.0, 0.0, 0.0), GeometryTypes.Vec(1.0, 1.0, 1.0))
c2 = GeometryTypes.HyperRectangle(GeometryTypes.Vec(3.0, 0.0, 0.0), GeometryTypes.Vec(1.0, 1.0, 1.0))

vis = Visualizer()
open(vis)
setobject!(vis[:c1], c1, MeshPhongMaterial(color=Colors.RGBA(1.0, 0.0, 0.0, 1.0)))
setobject!(vis, GeometryBasics.HyperRectangle(GeometryBasics.Vec(0., 0, 0), GeometryBasics.Vec(1., 1, 1)))
setobject!(vis, GeometryTypes.HyperRectangle(GeometryTypes.Vec(0., 0, 0), GeometryTypes.Vec(1., 1, 1)))

EnhancedGJK.dimension(::Type{GeometryBasics.AbstractSimplex{S, GeometryBasics.Vec{N, T}}}) where {N,T,S} = Val(N)
EnhancedGJK.dimension(::Type{M}) where {M <: GeometryBasics.AbstractMesh} = dimension(GeometryBasics.vertextype(M))
EnhancedGJK.dimension(::Type{GeometryBasics.Vec{N, T}}) where {N,T} = Val(N)
EnhancedGJK.dimension(::Type{GeometryBasics.Point{N, T}}) where {N,T} = Val(N)
EnhancedGJK.dimension(::Type{GeometryBasics.Simplex{M, T}}) where {M,T} = dimension(T)
EnhancedGJK.dimension(::Type{SVector{N, T}}) where {N,T} = Val(N)
# EnhancedGJK.dimension(::Type{GeometryBasics.FlexibleConvexHull{T}}) where {T} = dimension(T)
EnhancedGJK.dimension(::Type{GeometryBasics.HyperRectangle{N, T}}) where {N,T} = Val(N)
EnhancedGJK.dimension(::Type{GeometryBasics.Rect3D{T}}) where T = Val(3)

# EnhancedGJK.any_inside(c::GeometryBasics.Rect3D{T}) where T = GeometryBasics.origin(c)

function EnhancedGJK.any_inside(geometry)
    @show "hi"
    point = GeometryBasics.any_inside(geometry)
    Tagged(SVector(point))
end

EnhancedGJK.any_inside(c3)


function any_inside(mesh::GeometryBasics.AbstractMesh{GeometryBasics.Point{N, T}}) where {N,T}
    Tagged(SVector(first(GeometryBasics.vertices(mesh))))
end

function any_inside(point::GeometryBasics.Point)
    Tagged(SVector(point))
end


c3 = GeometryBasics.HyperRectangle(GeometryBasics.Vec(0.0, 0.0, 0.0), GeometryBasics.Vec(1.0, 1.0, 1.0))
c4 = GeometryBasics.HyperRectangle(GeometryBasics.Vec(3.0, 0.0, 0.0), GeometryBasics.Vec(1.0, 1.0, 1.0))

dimension(typeof(c3))

result1 = EnhancedGJK.gjk(c1, c2)
result2 = EnhancedGJK.gjk(c3, c4)

