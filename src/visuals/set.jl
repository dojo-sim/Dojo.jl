"""
    set_floor!(vis; x, y, z, origin, normal, color, tilepermeter, imagename, axis, grid)

    adds floor to visualization

    vis::Visualizer
    x: lateral position
    y: longitudinal position
    z: vertical position
	origin:: position of the center of the floor
    normal:: unit vector indicating the normal to the floor
    color: RGBA
    tilepermeter: scaling
    imagename: path to image
    axis: flag to turn on visualizer axis
    grid: flag to turn on visualizer grid
"""
function set_floor!(vis::Visualizer;
	    x=20.0,
	    y=20.0,
	    z=0.1,
	    origin=[0,0,0.0],
		normal=[0,0,1.0],
	    color=RGBA(0.5,0.5,0.5,1.0),
	    tilepermeter=1.0,
	    imagename="tile.png",
	    axis::Bool=false,
	    grid::Bool=false)
    image = PngImage(joinpath(module_dir(), "assets", imagename))
    repeat = Int.(ceil.(tilepermeter * [x, y]))
    texture = Texture(image=image, wrap=(1,1), repeat=(repeat[1],repeat[2]))
    mat = MeshPhongMaterial(map=texture)
    (color !== nothing) && (mat = MeshPhongMaterial(color=color))
    obj = HyperRectangle(Vec(-x/2, -y/2, -z), Vec(x, y, z))
    setobject!(vis[:floor], obj, mat)
	p = origin
	q = axes_pair_to_quaternion([0,0,1.], normal)
    settransform!(vis[:floor], MeshCat.compose(
		MeshCat.Translation(p...),
		MeshCat.LinearMap(q),
		))

    setvisible!(vis["/Axes"], axis)
    setvisible!(vis["/Grid"], grid)
    return nothing
end

"""
    set_surface!(vis; f, xlims, ylims, color, n)

    adds surface to visualization

    vis::Visualizer
    f: implicit function representing surface
    xlims: lateral domain for surface
    ylims: longitudinal domain for surface
    color: RGBA
    n: number of discretization points along each domain
"""
function set_surface!(vis::Visualizer, f::Any;
    xlims=[-20.0, 20.0],
    ylims=[-20.0, 20.0],
    color=RGBA(1.0, 1.0, 1.0, 0.0),
    n::Int=200)
    mesh = GeometryBasics.Mesh(f,
        HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2] - xlims[1], ylims[2] - ylims[1], 4.0)),
        Meshing.MarchingCubes(), samples=(n, n, Int(floor(n / 8))))
    setobject!(vis["surface"], mesh, MeshPhongMaterial(color=color))
    return nothing
end

"""
    set_light!(vis; ambient, fill, pointX, pointXshadow, direction)

    lighting conditions for visualization

    vis: Visualizer
    ambient: value for ambient lighting
    direction: positive or negative direction for light
"""
function set_light!(vis::Visualizer;
    ambient=0.35,
    direction::String="Positive")

    setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", ambient)
    setprop!(vis["/Lights/FillLight/<object>"], "intensity", 0.25)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "intensity", 0.85)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "castShadow", true)
    return nothing
end

"""
    set_camera!(vis; zoom, cam_pos)

    position and zoom for camera in visualization

    vis: Visualizer
    zoom: value for zoom
    cam_pos: position of camera
"""
function set_camera!(vis::Visualizer;
    zoom=1.0,
    cam_pos=nothing)

    camvis=vis["/Cameras/default/rotated/<object>"]
    setprop!(camvis, "zoom", zoom)
    (cam_pos !== nothing) && MeshCat.settransform!(camvis,
        MeshCat.compose(
            MeshCat.LinearMap(Dojo.RotX(-1/2 * pi)),
            MeshCat.Translation(cam_pos...),
        ))
    return nothing
end

"""
    set_arrow!(vis, origin, direction; color, shaft_radius, max_head_radius, scaling, name)

    adds an arrow object to scene 

    vis: Visualizer 
    origin: point defining arrow base 
    direction: vector defining arrow 
    color: RGBA 
    shaft_radius: dimension of arrow shaft 
    max_head_radius: dimension of arrow head base 
    scaling: parameter that scales the entire arrow 
    name: Symbol
"""
function set_arrow!(vis, origin, direction; 
    color=Colors.RGBA(1.0, 0.0, 0.0, 1.0), 
    shaft_radius=0.0125, 
    max_head_radius=0.025, 
    scaling=0.2, 
    name=:name)

    # create arrow
    force_vis = ArrowVisualizer(vis[name])
    setobject!(force_vis, MeshPhongMaterial(color=color))

    # direction 
    scaled_direction = scaling * direction 

    # set 
    settransform!(force_vis,
        Point(origin...),
        Vec(scaled_direction...),
        shaft_radius=shaft_radius,
        max_head_radius=max_head_radius)

end

function set_background!(vis::Visualizer; top_color=RGBA(1,1,1.0), bottom_color=RGBA(1,1,1.0))
    setprop!(vis["/Background"], "top_color", top_color)
    setprop!(vis["/Background"], "bottom_color", bottom_color)
end
