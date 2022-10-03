function main(directory, filename)
    read_filepath = joinpath(directory, filename)

    io = 0
    idx = 0
    for (i,line) in enumerate(eachline(read_filepath; keep=true))
        # detect new object
        if line[1] == 'o'
            (idx > 0) && close(io)
            idx += 1
            @show idx
            object_path = joinpath(directory, "convexes", "convex_$idx.obj")
            io = open(object_path, "w")
        end
        write(io, line)
    end
    close(io)
end


directory = joinpath(@__DIR__, "bluesoap", "decomposed")
filename = "bluesoap_level_3_decomposed.obj"

main(directory, filename)
