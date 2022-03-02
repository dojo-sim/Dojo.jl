# linear cone 
line_mat_lc = LineBasicMaterial(color=color_lc, linewidth=10.0)
points_lc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_lc.x[1])
    k = xt
    push!(points_lc, Point(k[1], k[2] + 0.04, k[3]))
end

setobject!(vis[:path_lc], Line(points_lc, line_mat_lc))

# nonlinear cone 
line_mat_nc = Dojo.LineBasicMaterial(color=color_nc, linewidth=25.0)
points_nc = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_nc.x[1])
    k = xt
    push!(points_nc, Point(k[1], k[2] - 0.04, k[3] + 0.0))
end
setobject!(vis[:path_nc], MeshCat.Line(points_nc, line_mat_nc))

# MuJoCo 
line_mat_mj = LineBasicMaterial(color=color_mj, linewidth=25.0)
points_mj = Vector{Point{3,Float64}}()
for (i, xt) in enumerate(storage_mj.x[1])
    k = xt
    push!(points_mj, Point(k[1], k[2], k[3]))
end
setobject!(vis[:path_mj], MeshCat.Line(points_mj, line_mat_mj))

settransform!(vis[:lc], MeshCat.Translation(0,+0.04,0))
settransform!(vis[:nc], MeshCat.Translation(0,-0.04,0))
settransform!(vis[:mj], MeshCat.Translation(0,+0.00,0))