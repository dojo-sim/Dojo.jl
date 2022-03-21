function normalize(x)
    mag = norm(x)
    if mag > 0.0 
        return x ./ mag 
    else 
        return ones(length(x)) ./ length(x)
    end
end

function ∂normalize∂x(x)
    mag = norm(x)
    n = length(x)
    if mag > 0.0 
        n = length(x)
        return 1.0 * I(n) ./ mag - x * transpose(x) ./ mag^3
    else 
        return 1.0 * I(n) # TODO: confirm this is good choice
    end
end

function ∂norm∂x(x)
    mag = norm(x) 
    if mag > 0.0 
        return x' ./ mag 
    else
        return ones(1, length(x))
    end
end