<img src="https://user-images.githubusercontent.com/4486578/59233116-77ec4a00-8c2a-11e9-87f3-012137752347.png" title="Solving a dual-valued linear system" align="center" width="100%"/>

# DualMatrixTools.jl

<p>
  <a href="https://doi.org/10.5281/zenodo.1493571">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1493571.svg" alt="DOI">
  </a>
  <a href="https://briochemc.github.io/DualMatrixTools.jl/stable">
    <img src=https://img.shields.io/badge/docs-stable-blue.svg>
  </a>
  <a href="https://travis-ci.com/briochemc/DualMatrixTools.jl">
    <img alt="Build Status" src="https://travis-ci.com/briochemc/DualMatrixTools.jl.svg?branch=master">
  </a>
  <a href="https://travis-ci.org/briochemc/DualMatrixTools.jl">
    <img alt="Build Status" src="https://travis-ci.org/briochemc/DualMatrixTools.jl.svg?branch=master">
  </a>
  <a href='https://coveralls.io/github/briochemc/DualMatrixTools.jl?branch=master'>
    <img src='https://coveralls.io/repos/github/briochemc/DualMatrixTools.jl/badge.svg?branch=master' alt='Coverage Status' />    
  </a>
  <a href="https://codecov.io/gh/briochemc/DualMatrixTools.jl">
    <img src="https://codecov.io/gh/briochemc/DualMatrixTools.jl/branch/master/graph/badge.svg" />
  </a>
  <a href="https://github.com/briochemc/DualMatrixTools.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
</p>

[![Linux](https://github.com/briochemc/DualMatrixTools.jl/workflows/Linux/badge.svg?branch=master)](https://github.com/briochemc/DualMatrixTools.jl/actions)
[![Mac OS X](https://github.com/briochemc/DualMatrixTools.jl/workflows/Mac%20OS%20X/badge.svg?branch=master)](https://github.com/briochemc/DualMatrixTools.jl/actions)
[![Windows](https://github.com/briochemc/DualMatrixTools.jl/workflows/Windows/badge.svg?branch=master)](https://github.com/briochemc/DualMatrixTools.jl/actions)

This package provides an overloaded `factorize` and `\` that work with dual-valued arrays.

It uses the dual type defined by the [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl) package.
The idea is that for a dual-valued matrix

*M* = *A* + *ϵ* *B*

its inverse is

*M*<sup>-1</sup> = (*I* - *ε* *A*<sup>-1</sup> *B*) *A*<sup>-1</sup>

Therefore, only the inverse of *A* is required to evaluate the inverse of *M*.
This package makes available a `DualFactors` type which containts (i) the factors of *A* and (ii) the non-real part, *B*.
It also overloads `factorize` to create an instance of `DualFactors` (when invoked with a dual-valued matrix), which can then be called with `\` to efficiently solve dual-valued linear systems of the type *M* *x* = *b*.

This package should be useful for autodifferentiation of functions that use `\`.
Note the same idea extends to hyper dual numbers (see the [HyperDualMatrixTools.jl](https://github.com/briochemc/HyperDualMatrixTools.jl) package).

## Usage

1. Create your dual-valued matrix `M`:
    ```julia
    julia> M = A + ε * B
    ```

2. Apply `\` to solve systems of the type `M * x = b`
    - without factorization:
        ```julia
        julia> x = M \ b
        ```

    - or better, with prior factorization:
        ```julia
        julia> Mf = factorize(M)

        julia> x = Mf \ b
        ```
        (This is better in case you want to solve for another `b`!)

## Advanced usage

In the context of iterative processes with multiple factorizations and forward and back substitutions, you may want to propagate dual-valued numbers while leveraging (potentially) the fact the real part of the matrices to be factorized remains the same throughout.
This package provides an in-place `factorize`, with a flag to update (or not) the factors.
Usage is straightforward.
By default, `factorize` does *not* update the factors
```julia
julia> factorize(Mf, M) # only Mf.B is updated
```
If you want to update the real-valued factors too, use
```julia
julia> factorize(Mf, M, update_factors=true) # Mf.B and Mf.Af are updated
```

## Citation

If you use this package, please cite it!
You can either directly use the bibtex entry, [CITATION.bib](CITATION.bib), or go to the [Zenodo record of the package](https://doi.org/10.5281/zenodo.1493571) and export the citation from there (the "Export" box at the bottom of that page).
