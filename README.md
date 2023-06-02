# ProbabilisticEchoInversion.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ElOceanografo.github.io/ProbabilisticEchoInversion.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ElOceanografo.github.io/ProbabilisticEchoInversion.jl/stable/)
[![Build Status](https://github.com/user/ProbabilisticEchoInversion.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/user/ProbabilisticEchoInversion.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for solving the "inverse problem" for multifrequency acoustic backscatter
in a Bayesian probabilistic framework. In other words, given observed echoes from
one or more types of scatterers, this package will help you infer 1) **what** they were, 
2) **how many** of them were present, and 3) **how sure** you can be about (1) and (2). We
call this approach Automatic Probabilistic Echo Solving, or APES.

Refer to the [documentation](https://ElOceanografo.github.io/ProbabilisticEchoInversion.jl/stable/) 
for more details on how to use the package. Several fully-worked examples on both simulated and 
real fisheries acoustic data can also be found in the 
[APES Examples](https://github.com/ElOceanografo/APESExamples) repository.

If you use this package in your own research, please cite the following publication:

Urmy, De Robertis, and Bassett (2023). A Bayesian inverse approach to identify and quantify
organisms from fisheries acoustic data. ICES Journal of Marine Science x (x), doi:x
