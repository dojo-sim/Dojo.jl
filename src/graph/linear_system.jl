"""
    full_matrix(system) 

    returns matrix from a simulation step's linear system (i.e., A from "Ax = b")

    system: System
"""
full_matrix(system::System) = full_matrix(system.matrix_entries, system.dimrow, system.dimcol)

function full_matrix(matrix_entries::SparseMatrixCSC, dimrow, dimcol)
    range_row = [1:dimrow[1]]
    for (i,dim) in enumerate(collect(Iterators.rest(dimrow, 2)))
        push!(range_row, sum(dimrow[1:i])+1:sum(dimrow[1:i])+dim)
    end
    range_col = [1:dimcol[1]]
    for (i,dim) in enumerate(collect(Iterators.rest(dimcol, 2)))
        push!(range_col, sum(dimcol[1:i])+1:sum(dimcol[1:i])+dim)
    end

    A = zeros(sum(dimrow),sum(dimcol))

    for (i,row) in enumerate(matrix_entries.rowval)
        col = findfirst(x->i<x,matrix_entries.colptr)-1
            A[range_row[row],range_col[col]] = matrix_entries[row,col].value
    end
    return A
end

# assemble gradient (i.e., b from "Ax = b")
"""
    full_vector(system) 

    returns residual vector from a simulation step's linear system (i.e., b from "Ax = b")

    system: System
"""
full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
