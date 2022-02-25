function module_dir()
    return joinpath(@__DIR__, "..", "..")
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
