using Test

using DualMatrixTools
using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

@testset "Testing DualMatrixTools" begin
    # Chose a size for matrices
    n = 10

    @testset "Testing full matrices" begin
        # Create a real-valued random vector and matrix
        y = randn(n)
        A = randn(n, n)

        # Create a hyperdual-valued random vector and matrix
        x = randn(n, 2) * [1.0, ε]
        B = randn(n, n)
        M = A + ε * B

        # Check that `\` works without factorization
        @test A \ (A * x) ≈ x
        @test M \ (M * x) ≈ x
        @test A \ (A * y) ≈ y
        @test M \ (M * y) ≈ y
        @test A * (A \ x) ≈ x
        @test M * (M \ x) ≈ x
        @test A * (A \ y) ≈ y
        @test M * (M \ y) ≈ y

        # Check that `\` works with factorization
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
    end

    @testset "Testing sparse matrices" begin
        # Create a real-valued sparse matrix
        A = sparse(randn(n, n))

        # Create a hyperdual-valued sparse matrix
        B = sparse(randn(n, n))
        M = A + ε * B

        # Check that `\` works without factorization
        @test A \ (A * x) ≈ x
        @test M \ (M * x) ≈ x
        @test A \ (A * y) ≈ y
        @test M \ (M * y) ≈ y
        @test A * (A \ x) ≈ x
        @test M * (M \ x) ≈ x
        @test A * (A \ y) ≈ y
        @test M * (M \ y) ≈ y

        # Check that `\` works with factorization
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
    end

    @testset "Testing isapprox function" begin
        A2 = A .+ 0ε
        @test A2 ≈ A
        @test A ≈ A2
        @test (1.0 + 0ε) ≈ 1.0
        @test 1.0 ≈ (1.0 + 0ε)
    end
end
