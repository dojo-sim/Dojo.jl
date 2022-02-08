using FFMPEG
using Meshing

function set_floor!(vis::Visualizer; x=20.0, y=20.0, z=0.1, alt=0.0, color=RGBA(0.5,0.5,0.5,1.0),
        tilepermeter=1.0, imagename="tile.png")
    image = PngImage(joinpath(module_dir(), "assets", imagename))
    repeat = Int.(ceil.(tilepermeter * [x,y]))
    texture = Texture(image=image, wrap=(1,1), repeat=(repeat[1],repeat[2]))
    mat = MeshPhongMaterial(map=texture)
    (color != nothing) && (mat = MeshPhongMaterial(color=color))
    obj = HyperRectangle(Vec(0., 0, 0), Vec(x, y, z))
    setobject!(vis[:floor], obj, mat)
    settransform!(vis[:floor], MeshCat.Translation(-x/2, -y/2, -z+alt))
    return nothing
end

function set_surface!(vis::Visualizer, f::Any; xlims = [-20.0, 20.0],
    ylims = [-20.0, 20.0], color=RGBA(1.0, 1.0, 1.0, 0.0), n::Int=200)
    mesh = GeometryBasics.Mesh(f,
        HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2]-xlims[1], ylims[2]-ylims[1], 4.0)),
        Meshing.MarchingCubes(), samples=(n, n, Int(floor(n/8))))
    setobject!(vis["surface"], mesh, MeshPhongMaterial(color=color))
    return nothing
end

function set_light!(vis::Visualizer; ambient=0.35, fill=0.25, pointX=0.85,
        pointXshadow::Bool=true, direction::String="Positive")
    setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", ambient)
    setprop!(vis["/Lights/FillLight/<object>"], "intensity", 0.25)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "intensity", 0.85)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "castShadow", true)
    return nothing
end

"""
    The camera always point towards the origin of the frame, you can choose its
    position with `cam_pos` and the `zoom`.
"""
function set_camera!(vis::Visualizer; zoom=1.0, cam_pos=nothing)
    camvis = vis["/Cameras/default/rotated/<object>"]
    setprop!(camvis, "zoom", zoom)
    (cam_pos != nothing) && MeshCat.settransform!(camvis,
        MeshCat.compose(
            MeshCat.LinearMap(Rotations.RotX(-1/2 * pi)),
            MeshCat.Translation(cam_pos...),
        ))
    return nothing
end

function convert_video_to_gif(video_file_path::String, output_path::String="output.gif";
    framerate::Int=30, start_time=0., duration=1e3, overwrite=false, width::Int=1080, height::Int=-2, hq_colors::Bool=false)
    output_path = abspath(output_path)

    if !isfile(video_file_path)
        error("Could not find the input file $video_file_path")
    end
    if isfile(output_path) && !overwrite
        error("The output path $output_path already exists. To overwrite that file, you can pass `overwrite=true` to this function")
    end

    mktempdir() do tmpdir
        # run(MeshCat.unpack_cmd(video_file_path, tmpdir, ".mp4", nothing)) # unpack the .tar file
        # cmd = ["-r", string(framerate), "-i", "%07d.png", "-vcodec", "libx264", "-preset", "slow", "-crf", "18"]
        color_map = hq_colors ?
            "[0:v] fps=$framerate, scale=$width:$height,split [a][b];[a] palettegen=stats_mode=single [p];[b][p] paletteuse=new=1" :
            "[0:v] fps=$framerate, scale=$width:$height,split [a][b];[a] palettegen [p];[b][p] paletteuse"
        cmd = ["-ss", string(start_time), "-t", string(duration), "-i", video_file_path, "-filter_complex", color_map]
        if overwrite
            push!(cmd, "-y")
        end
        push!(cmd, output_path)

        cd(tmpdir) do
            FFMPEG.exe(cmd...)
        end
    end
    @info("Saved output as $output_path")
    return output_path
end

function convert_frames_to_video_and_gif(filename, overwrite::Bool=true)
    MeshCat.convert_frames_to_video(
        homedir() * "/Downloads/$filename.tar",
        homedir() * "/Documents/video/$filename.mp4", overwrite=overwrite)

    convert_video_to_gif(
        homedir() * "/Documents/video/$filename.mp4",
        homedir() * "/Documents/video/$filename.gif", overwrite=overwrite)
    return nothing
end
