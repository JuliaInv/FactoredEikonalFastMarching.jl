# FactoredEikonalFastMarching.jl
Julia Package for solving the factored eikonal equation of a regular rectangular mesh using the fast marching algorithm.

# Requirements

This package is intended to use with julia versions 0.4.x.

This package is an add-on for jInv, which needs to be installed. This is for having a Mesh module.

# Installation


Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/FactoredEikonalFastMarching.jl","FactoredEikonalFastMarching")
Pkg.test("DivSigGrad")
