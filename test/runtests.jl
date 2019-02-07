using Test

using DualMatrixTools
using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

# Check definition of dual constant
@test ε == dual(0.0, 1.0)

# Chose a size for matrices
n = 10

# Create a real-valued random vector and matrix
y = randn(n)
A = randn(n, n)

# Create a hyperdual-valued random vector and matrix
x = randn(n, 2) * [1.0, ε]
B = randn(n, n)
M = A + ε * B

# Check that `\` works for full matrices
Af = factorize(A)
Mf = factorize(M)
@test Af \ (A * x) ≈ x
@test Mf \ (M * x) ≈ x
@test Af \ (A * y) ≈ y
@test Mf \ (M * y) ≈ y
@test A * (Af \ x) ≈ x
@test M * (Mf \ x) ≈ x
@test A * (Af \ y) ≈ y
@test M * (Mf \ y) ≈ y

# Create a real-valued sparse matrix
A = sparse(randn(n, n))

# Create a hyperdual-valued sparse matrix
B = sparse(randn(n, n))
M = A + ε * B

# Check that `\` works for sparse matrices
Af = factorize(A)
Mf = factorize(M)
@test Af \ (A * x) ≈ x
@test Mf \ (M * x) ≈ x
@test Af \ (A * y) ≈ y
@test Mf \ (M * y) ≈ y
@test A * (Af \ x) ≈ x
@test M * (Mf \ x) ≈ x
@test A * (Af \ y) ≈ y
@test M * (Mf \ y) ≈ y

# Check that isapprox is used in all the ways (for code coverage)
A2 = A + ε * B * 0
@test A2 ≈ A
@test A ≈ A2
@test (1.0 + 0ε) ≈ 1.0
@test 1.0 ≈ (1.0 + 0ε)
