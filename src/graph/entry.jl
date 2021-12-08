mutable struct Entry{ET}
    value::ET
    isinverted::Bool

    function Entry{T}(dims...; static = true) where T
        static ? value = szeros(T,dims...) : value = zeros(T,dims...)
        new{typeof(value)}(value, false)
    end

    Entry(dims...; static = true) = Entry{Float64}(dims...; static = static)
end

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, entry::Entry)
#     println(io, "Entry with value:")
#     show(io, mime, entry.value)
# end

function Base.zero(::Entry{ET}) where ET
    dims = [ET.parameters[1].parameters...]
    return Entry{ET.parameters[2]}(dims...)
end
function Base.zero(::Type{Entry{ET}}) where ET
    dims = [ET.parameters[1].parameters...]
    return Entry{ET.parameters[2]}(dims...)
end

function randomize!(entry::Entry)
    value = entry.value
    entry.value = randn(eltype(value), size(value))
end
