# FactoredEikonalFastMarching.jl
Julia Package for solving the factored eikonal equation on a regular rectangular mesh using the fast marching algorithm.

Based on the following paper (please cite if you are using the package):

Eran Treister and Eldad Haber, A fast marching algorithm for the factored eikonal equation, Under review.

# Requirements

This package is intended to use with julia versions 0.4.x.

This package is an add-on for jInv, which needs to be installed. This is for having a Mesh module.

# Installation

```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/FactoredEikonalFastMarching.jl","FactoredEikonalFastMarching")
Pkg.test("FactoredEikonalFastMarching")
```

# Examples

Under "examples/runExperiments.jl" you can find the experiments that were shown in the paper above. 


