using MeshCat
using Colors

using EnhancedGJK 
using CoordinateTransformations
using StaticArrays

using GeometryTypes
using GeometryBasics

r1 = GeometryTypes.HyperRectangle(GeometryTypes.Vec(0.0, 0.0, 0.0), GeometryTypes.Vec(1.0, 1.0, 1.0))
r2 = GeometryTypes.HyperRectangle(GeometryTypes.Vec(0.0, 0.0, 0.0), GeometryTypes.Vec(1.0, 1.0, 1.0))

r3 = GeometryBasics.HyperRectangle(GeometryBasics.Vec(0.0, 0.0, 0.0), GeometryBasics.Vec(1.0, 1.0, 1.0))
r4 = GeometryBasics.HyperRectangle(GeometryBasics.Vec(0.0, 0.0, 0.0), GeometryBasics.Vec(1.0, 1.0, 1.0))

s1 = GeometryTypes.HyperSphere(GeometryTypes.Point(0.0, 0.0, 0.0, 0.0), 1.0)
s2 = GeometryBasics.HyperSphere(GeometryBasics.Point(0.0, 0.0, 0.0, 0.0), 1.0)


vis = Visualizer()
open(vis)


result1 = EnhancedGJK.gjk(r1, r2)
result2 = GeometryTypes.gjk(s1, s2)

