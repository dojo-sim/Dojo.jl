

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

function plot_cone(s, Δs; plt = plot(), xlabel = "", linewidth = 3.0, show = true, zoom = false,
		xlims = (-1.3,1.3), ylims = (-1.3,1.3))
	# Current point
	px = 1/s[1]*s[2]
	py = 1/s[1]*s[3]
	window = 1.3 * max(1-norm(s[2:3]/s[1]), norm(Δs))
	if zoom
		xlims = (-window + px, +window + px)
		ylims = (-window + py, +window + py)
	end
	θ = 0:0.001:2π
	Xc = cos.(θ)
	Yc = sin.(θ)
	plot!(plt, Xc, Yc, legend = false, aspectratio = 1.0, xlims = xlims, ylims = ylims, xlabel = xlabel, xtick = [], ytick = [])
	scatter!(plt, [px], [py], markersize = 5.0)
	α = 0:0.05:1.0
	Sα = [s + αi*Δs for αi in α]
	Xα = [1/s[1]*s[2] for s in Sα]
	Yα = [1/s[1]*s[3] for s in Sα]
	plot!(plt, Xα, Yα, linewidth = linewidth)
	show && display(plt)
	return plt
end
