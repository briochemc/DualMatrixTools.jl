using Test

using DualMatrixTools
using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

@testset "Testing DualMatrixTools" begin
    # Chose a size for matrices
    n = 10
    y = randn(n)               # real-valued vector
    x = randn(n, 2) * [1.0, ε] # dual-valued vector
    A = randn(n, n)            # real-valued matrix
    B = randn(n, n)            # real-valued matrix

    @testset "Testing full matrices" begin
        # Create a dual-valued matrix
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
        A = sparse(A)

        # Create a hyperdual-valued sparse matrix
        B = sparse(B)
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

    @testset "Testing adjoint and transpose" for f in (:adjoint, :transpose)
        # Create a real-valued sparse matrix
        A = sparse(A)

        # Create a hyperdual-valued sparse matrix
        B = sparse(B)
        M = A + ε * B

        # Adjoints / transposes
        y2 = randn(n)               # real-valued vector
        x2 = randn(n, 2) * [1.0, ε] # dual-valued vector

        @eval( :( ay  = $f(y)  ) )
        @eval( :( ax  = $f(x)  ) )
        @eval( :( ay2 = $f(y2) ) )
        @eval( :( ax2 = $f(x2) ) )
        @eval( :( aA  = $f(A)  ) )
        @eval( :( aM  = $f(M)  ) )

        # Check that `\` works with adjoints
        @test (ax2 * (A \ x)) ≈ (x * (aA \ ax2))
        @test (ay2 * (A \ x)) ≈ (x * (aA \ ay2))
        @test (ax2 * (M \ x)) ≈ (x * (aM \ ax2))
        @test (ay2 * (M \ x)) ≈ (x * (aM \ ay2))
        @test (ax2 * (A \ y)) ≈ (y * (aA \ ax2))
        @test (ay2 * (A \ y)) ≈ (y * (aA \ ay2))
        @test (ax2 * (M \ y)) ≈ (y * (aM \ ax2))
        @test (ay2 * (M \ y)) ≈ (y * (aM \ ay2))

        # Check that `\` works with adjoints of factorized versions
        Af = factorize(A)
        Mf = factorize(M)
        @eval( :( aAf  = $f(Af)  ) )
        @eval( :( aMf  = $f(Mf)  ) )
        @test (ax2 * (Af \ x)) ≈ (x * (aAf \ ax2))
        @test (ay2 * (Af \ x)) ≈ (x * (aAf \ ay2))
        @test (ax2 * (Mf \ x)) ≈ (x * (aMf \ ax2))
        @test (ay2 * (Mf \ x)) ≈ (x * (aMf \ ay2))
        @test (ax2 * (Af \ y)) ≈ (y * (aAf \ ax2))
        @test (ay2 * (Af \ y)) ≈ (y * (aAf \ ay2))
        @test (ax2 * (Mf \ y)) ≈ (y * (aMf \ ax2))
        @test (ay2 * (Mf \ y)) ≈ (y * (aMf \ ay2))
    end

    @testset "Testing isapprox function" begin
        A2 = A .+ 0ε
        @test A2 ≈ A
        @test A ≈ A2
        @test (1.0 + 0ε) ≈ 1.0
        @test 1.0 ≈ (1.0 + 0ε)
    end
end
