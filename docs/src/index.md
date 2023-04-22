```@meta
CurrentModule = ProbabilisticEchoInversion
```

# ProbabilisticEchoInversion

Documentation for [ProbabilisticEchoInversion](https://github.com/user/ProbabilisticEchoInversion.jl).

## Introduction

Welcome to the documentation for ProbabilisticEchoInversion.jl!

This package is designed to solve the "inverse problem" for acoustic backscatter at
multiple frequencies in a Bayesian framework. In other words, given observed echoes from
one or more types of scatterers, this package will help you infer 1) **what** they were, 
2) **how many** of them were present, and 3) **how sure** you can be about (1) and (2). We
call this approach Automatic Probabilistic Echo Solving, or APES.

In a nutshell, you arrange your acoustic data in a multidimensional `DimArray` from 
DimensionalData.jl.[^1] For a typical downard-looking echosounder on a moving ship, your 
data would have three dimensions, depth x distance x frequency.

You then specify the inverse problem you want to solve using the modeling language Turing.jl, 
and then use ProbabilisticEchoInversion to apply it to each depth/distance cell in your data


[^1] (`DimArray`s behave like normal Julia arrays, but also 
let you index with named dimensions and have a bunch of nice functionality for subsetting
and plotting.) 

## Tutorial

## API

```@index
```

```@autodocs
Modules = [ProbabilisticEchoInversion]
```
