"""
    Entry{ET}

    value: data 
    isinverted: flag indicating if matrix value has been inverted
"""
mutable struct Entry{ET}
    value::ET
    isinverted::Bool

    function Entry{T}(dims...; 
        static = true) where T
        static ? value = szeros(T,dims...) : value = zeros(T,dims...)
        new{typeof(value)}(value, false)
    end

    Entry(dims...; static = true) = Entry{Float64}(dims...; static = static)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, entry::Entry)
    summary(io, entry)
    println(io, "")
    println(io, "value:      "*string(entry.value))
    println(io, "isinverted: "*string(entry.isinverted))
end

function Base.zero(::Entry{ET}) where ET
    dims = [ET.parameters[1].parameters...]
    return Entry{ET.parameters[2]}(dims...)
end

Base.zero(::Type{Dojo.Entry}) = nothing
