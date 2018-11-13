module DualMatrixTools

using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

import LinearAlgebra.factorize
import Base.\
import Base.isapprox

"""
    DualFactors

Container type to work efficiently with backslash on dual-valued sparse matrices.

`factorize(M)` will create an instance containing
- `Af = factorize(realpart.(M))` — the factors of the real part
- `B = dualpart.(M)` — the dual part
for a dual-valued matrix `M`.

This is because only the factors of the real part are needed when solving a linear system of the type ``M x = b`` for a dual-valued matrix ``M = A + \\varepsilon B``.
In fact, the inverse of ``M`` is given by
``M^{-1} = (I - \\varepsilon A^{-1} B) A^{-1}``.
"""
mutable struct DualFactors
    Af::Factorization # the factors of the real part
    B                 # the ε₁ part
end


"""
    factorize(M::Array{Dual128,2})

Efficient factorization of dual-valued matrices.
See `DualFactors` for details.
"""
function factorize(M::Array{Dual128,2})
    return DualFactors(factorize(realpart.(M)), dualpart.(M))
end

"""
    factorize(M::SparseMatrixCSC{Dual128,Int64})

Efficient factorization of dual-valued sparse matrices.
See `DualFactors` for details.
"""
function factorize(M::SparseMatrixCSC{Dual128,Int64})
    return DualFactors(factorize(realpart.(M)), dualpart.(M))
end

"""
    \\(M::DualFactors, y::AbstractVecOrMat{Float64})

Backsubstitution for `DualFactors`.
See `DualFactors` for details.
"""
function \(M::DualFactors, y::AbstractVecOrMat{Float64})
    A, B = M.Af, M.B
    A⁻¹y = A \ y
    A⁻¹BA⁻¹y = A \ (B * A⁻¹y)
    return A⁻¹y - ε * A⁻¹BA⁻¹y
end

"""
    \\(M::DualFactors, y::AbstractVecOrMat{Dual128})

Backsubstitution for `DualFactors`.
See `DualFactors` for details.
"""
function \(M::DualFactors, y::AbstractVecOrMat{Dual128})
    a, b = realpart.(y), dualpart.(y)
    A, B = M.Af, M.B
    A⁻¹a = A \ a
    A⁻¹BA⁻¹a = A \ (B * A⁻¹a)
    A⁻¹b = A \ b
    return A⁻¹a - ε * A⁻¹BA⁻¹a + ε * A⁻¹b
end

"""
    \\(Af::Factorization{Float64}, y::AbstractVecOrMat{Dual128})

Backsubstitution for Dual-valued RHS.
"""
function \(Af::Factorization{Float64}, y::AbstractVecOrMat{Dual128})
    return (Af \ realpart.(y)) + ε * (Af \ dualpart.(y))
end

function isapprox(x::AbstractVecOrMat{Dual128}, y::AbstractVecOrMat{Dual128})
    bigx = [realpart.(x) dualpart.(x)]
    bigy = [realpart.(y) dualpart.(y)]
    return isapprox(bigx, bigy)
end
isapprox(x::AbstractVecOrMat, y::AbstractVecOrMat{Dual128}) = isapprox(dual.(x), y)
isapprox(x::AbstractVecOrMat{Dual128}, y::AbstractVecOrMat) = isapprox(x, dual.(y))
function isapprox(x::Dual128, y::Dual128)
    bigx = [realpart(x) dualpart(x)]
    bigy = [realpart(y) dualpart(y)]
    return isapprox(bigx, bigy)
end
isapprox(x::Float64, y::Dual128) = isapprox(dual(x), y)
isapprox(x::Dual128, y::Float64) = isapprox(x, dual(y))

export DualFactors, factorize, \

end # module
