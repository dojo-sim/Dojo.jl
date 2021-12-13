function deleteat(M::Array, i1::Integer, i2::Integer)
    return [M[1:i1 - 1,1:i2 - 1] M[1:i1 - 1,i2 + 1:end];M[i1 + 1:end,1:i2 - 1] M[i1 + 1:end,i2 + 1:end]]
end

deleteat(M::Array,i::Integer) = deleteat(M, i, i)

function orthogonalrows(axis::AbstractVector)
    if norm(axis) > 0
        axis = normalize(axis)
    end
    A = svd(skew(axis)).Vt
    inds = SA[1; 2; 3]
    V1 = A[1,inds]'
    V2 = A[2,inds]'
    V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

    return V1, V2, V3
end

function orthogonalcols(axis::AbstractVector)
    V1, V2, V3 = orthogonalrows(axis)
    return V1', V2', V3'
end

@inline function skewplusdiag(v::AbstractVector{T},w::T) where T
    SA[
         w    -v[3]  v[2]
         v[3]  w    -v[1]
        -v[2]  v[1]  w
    ]
end

function getfieldnumber(obj)
    i = 1
    while true
        !isdefined(obj, i) ? break : (i+=1)
    end
    return i-1
end

@inline offsetrange(offset, length) = (offset-1)*length+1:offset*length
@inline offsetrange(offset, length, totallength, inneroffset) = (offset-1)*totallength+(inneroffset-1)*length+1:(offset-1)*totallength+inneroffset*length



function scn(a::Number; digits::Int=1, exp_digits::Int=1)
	(typeof(a) <: Float64) ? nothing : return nothing
end

function scn(a::Float64; digits::Int=1, exp_digits::Int=1)
	isnan(a) && return " NaN" * " "^(digits + exp_digits)
	@assert digits >= 0
    # a = m x 10^e
    if a == 0
        e = 0
        m = 0.0
    elseif a == Inf
		return " Inf"
	elseif a == -Inf
		return "-Inf"
	else
        e = Int(floor(log(abs(a))/log(10)))
        m = a*exp(-e*log(10))
    end

    m = round(m, digits=digits)
	if m == 10.0
		m = 1.0
		e += 1
	end
    if digits == 0
        m = Int(floor(m))
		strm = string(m)
	else
		strm = string(m)
		is_neg = m < 0.
		strm = strm*"0"^max(0, 2+digits+is_neg-length(strm))
    end
    sgn = a >= 0 ? " " : ""
    sgne = e >= 0 ? "+" : "-"

	stre = string(abs(e))
	stre = "0"^max(0, exp_digits - length(stre)) * stre
    return "$sgn$(strm)e$sgne$(stre)"
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

function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# function Dojo.mean(x)
# 	return sum(x) / length(x)
# end

function fdjac(f, x; δ = 1e-5)
    n = length(f(x))
    m = length(x)
    jac = zeros(n, m)
    for i = 1:m
        xp = deepcopy(x)
        xm = deepcopy(x)
        xp[i] += δ
        xm[i] -= δ
        jac[:,i] = (f(xp) - f(xm)) / (2δ)
    end
    return jac
end

# useful for python visualizer
wait_for_server = MeshCat.wait_for_server
