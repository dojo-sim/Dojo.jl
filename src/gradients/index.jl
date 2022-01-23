
################################################################################
# Index and Dimensions
################################################################################

function Fz_indices(Nb::Int)
    return vcat([13*(i-1) .+ [4:6; 11:13] for i = 1:Nb]...)
end

function datadim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; attjac::Bool = true) where {T,Nn,Ne,Nb,Ni}
    d = 0
    d += 12Nb
    !attjac && (d += Nb)
    d += controldim(mechanism)
    return d
end

function soldim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    d = 0
    d += 6Nb
    d += eqcdim(mechanism)
    d += ineqcdim(mechanism)
    return d
end

function controldim(eqc::EqualityConstraint{T,N,Nc,Cs}; ignore_floating_base::Bool = false) where {T,N,Nc,Cs}
    ignore_floating_base && (N == 0) && return 0

    N̄ = 0
    for (i, joint) in enumerate(eqc.constraints)
        N̄ += controldim(joint)
    end

    return N̄
end

function controldim(joint::Joint{T,N}) where {T,N}
    return 3 - N
end

function controldim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; ignore_floating_base::Bool = false) where {T,Nn,Ne,Nb,Ni}
    nu = 0
    for eqc in mechanism.eqconstraints
        nu += controldim(eqc, ignore_floating_base = ignore_floating_base)
    end
    return nu
end

function minCoordDim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nx = 0
    free_rot_base = false # we are going to check if the link attached to the base has free orientation
    nx = 2 * controldim(mechanism, ignore_floating_base = false)
    free_rot_base && (nx += 1)
    return nx
end

maxCoordDim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni} = 13Nb

function ineqcdim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nineqcs = 0
    for ineqc in mechanism.ineqconstraints
        nineqcs += length(ineqc)
    end
    return nineqcs
end

function eqcdim(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    neqcs = 0
    for eqc in mechanism.eqconstraints
        neqcs += length(eqc)
    end
    return neqcs
end
