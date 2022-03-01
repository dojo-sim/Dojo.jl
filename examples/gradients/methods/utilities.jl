function write_csv(names, data, filename)
    open(joinpath(@__DIR__, filename), "w") do f
        write(f, join(names, ",") * "\n")
        for tup in data
            write(f, join(tup, ",") * "\n")
        end
    end
    return nothing
end

function normal_sample(μ, Σ)
    n = length(μ)
    x = μ + sqrt(Σ) * randn(n)
    return x
end