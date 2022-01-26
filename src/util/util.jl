function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

function deleteat(M::Array, i1::Integer, i2::Integer)
    return [M[1:i1 - 1,1:i2 - 1] M[1:i1 - 1,i2 + 1:end];M[i1 + 1:end,1:i2 - 1] M[i1 + 1:end,i2 + 1:end]]
end

deleteat(M::Array,i::Integer) = deleteat(M, i, i)

function orthogonal_rows(axis::AbstractVector)
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
    V1, V2, V3 = orthogonal_rows(axis)
    return V1', V2', V3'
end

function get_fieldnumber(obj)
    i = 1
    while true
        !isdefined(obj, i) ? break : (i+=1)
    end
    return i-1
end

@inline offset_range(offset, length) = (offset-1)*length+1:offset*length
@inline offset_range(offset, length, totallength, inneroffset) = (offset-1)*totallength+(inneroffset-1)*length+1:(offset-1)*totallength+inneroffset*length

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

# useful for python visualizer
wait_for_server = MeshCat.wait_for_server
