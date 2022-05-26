# LDU factorization for unsymmetric systems
function ldu_factorization_acyclic!(diagonal_v, offdiagonal_l, diagonal_c, offdiagonal_u, diagonal_inverse_c)
    if diagonal_inverse_c.isinverted
        invdiagonal_c = diagonal_inverse_c.value
    else
        invdiagonal_c = inv(diagonal_c.value)
        diagonal_inverse_c.value = invdiagonal_c
        diagonal_inverse_c.isinverted = true
    end
    offdiagonal_l.value = offdiagonal_l.value * invdiagonal_c
    offdiagonal_u.value = invdiagonal_c * offdiagonal_u.value

    if length(diagonal_c.value) > 0 # mutliplication of matrices of size 0 is ambiguous; e.g. 6x0 * 0x0 * 0x6 = ???
        diagonal_v.value -= offdiagonal_l.value * diagonal_c.value * offdiagonal_u.value
    end
    return
end

function ldu_factorization_cyclic!(entry_lu, offdiagonal_lu, diagonal_c, offdiagonal_ul)
    entry_lu.value -= offdiagonal_lu.value * diagonal_c.value * offdiagonal_ul.value
    return
end

function ldu_factorization!(system)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children

    for v in system.dfs_list
        for c in cyclic_children[v]
            for cc in cyclic_children[v]
                cc == c && break
                (cc ∉ children(system,c) && cc ∉ cyclic_children[c]) && continue
                ldu_factorization_cyclic!(matrix_entries[v, c], matrix_entries[v, cc], matrix_entries[cc, cc], matrix_entries[cc, c])
                ldu_factorization_cyclic!(matrix_entries[c, v], matrix_entries[c, cc], matrix_entries[cc, cc], matrix_entries[cc, v])
            end
            ldu_factorization_acyclic!(matrix_entries[v, v], matrix_entries[v, c], matrix_entries[c, c], matrix_entries[c, v], diagonal_inverses[c])
        end
        for c in acyclic_children[v]
            ldu_factorization_acyclic!(matrix_entries[v, v], matrix_entries[v, c], matrix_entries[c, c], matrix_entries[c, v], diagonal_inverses[c])
        end
    end

    return
end

function ldu_backsubstitution_l!(vector_v, offdiagonal, vector_c)
    vector_v.value -= offdiagonal.value * vector_c.value
    return
end

function ldu_backsubstitution_u!(vector_v, offdiagonal, vector_p)
    vector_v.value -= offdiagonal.value * vector_p.value
    return
end

function ldu_backsubstitution_d!(vector, diagonal, diagonal_inverse)
    if diagonal_inverse.isinverted
        vector.value = diagonal_inverse.value * vector.value
    else
        vector.value = diagonal.value \ vector.value
    end
    diagonal_inverse.isinverted = false
    return
end

function ldu_backsubstitution!(system)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    vector_entries = system.vector_entries
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children
    parents = system.parents
    dfs_list = system.dfs_list

    for v in dfs_list
        for c in cyclic_children[v]
            ldu_backsubstitution_l!(vector_entries[v], matrix_entries[v, c], vector_entries[c])
        end
        for c in acyclic_children[v]
            ldu_backsubstitution_l!(vector_entries[v], matrix_entries[v, c], vector_entries[c])
        end
    end

    for v in reverse(dfs_list)
        ldu_backsubstitution_d!(vector_entries[v], matrix_entries[v, v], diagonal_inverses[v])
        for p in parents[v]
            ldu_backsubstitution_u!(vector_entries[v], matrix_entries[v, p], vector_entries[p])
        end
    end
end
