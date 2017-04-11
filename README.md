[![Build Status](https://travis-ci.org/JuliaInv/FactoredEikonalFastMarching.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/FactoredEikonalFastMarching.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/FactoredEikonalFastMarching.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/FactoredEikonalFastMarching.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/9pqt8ragr0icc9ss?svg=true)](https://ci.appveyor.com/project/lruthotto/factoredeikonalfastmarching-jl)


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

Under "examples/simpleExample.jl" you can find how to set a simple travel time calculation.

Under "examples/runExperiments.jl" you can find the experiments that were shown in the paper above. 


