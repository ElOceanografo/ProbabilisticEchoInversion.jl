```@meta
CurrentModule = ProbabilisticEchoInversion
```

# ProbabilisticEchoInversion.jl
## Introduction

Welcome to the documentation for ocumentation for 
[ProbabilisticEchoInversion](https://github.com/user/ProbabilisticEchoInversion.jl)!

This package is designed to solve the "inverse problem" for acoustic backscatter at
multiple frequencies in a Bayesian statistical framework. In other words, given observed 
echoes from one or more types of scatterers, this package will help you infer

1) **What** they were,
2) **How many** of them were present, and
3) **How sure** you can be about (1) and (2). 

We call this approach Automatic Probabilistic Echo Solving, or APES.

This documentation provides a short introduction to the problem and the general
approach, as well as a simple tutorial on how to use the package. It assumes
basic familiarity with the principles of fisheries acoustics and Bayesian
statistical modeling (e.g. experience with Turing.jl, Stan, or JAGS). While the package
is written in [Julia](https://julialang.org/), we don't assume knowledge of that language, and you should
hopefully be able to follow along as long as you are familiar with scientific programming
in a similar scripting language like R, Python, or Matlab. Extensive resources for learning
Julia are available at the [official website](https://julialang.org/learning/).

### Background
The typical approach in fisheries acoustics is to use relative differences in backscatter
at multiple frequencies to classify regions of the water column as one thing or another
(fish vs. zooplankton, large vs. small fish, etc.). These frequency responses are usually
calculated by taking the difference between frequencies in the decibel domain (dB 
differencing). Once the echogram has been classified using multiple frequencies,
scatterer density is estimated by echo-integrating at a single frequency using the 
relationship 

$s_v = \langle \sigma_{bs} \rangle n,$ (1)

where $s_v$ is the volume backscattering coefficient, $\langle \sigma_{bs} \rangle$ is the
average backscattering cross section of one scatterer, and $n$ is their numerical density.

### The Inverse Approach
The inverse approach relies on the same theory and assumptions, but instead of classifying 
first using multiple frequencies and then integrating using only one, it does both at 
the 
same time using all available frequencies. Mathematically, if we have a vector $\mathbf{s}_v$
of backscatter at multiple frequencies and a matrix $\Sigma$ of backscattering cross-sections,
where the $i,j^{th}$ entry holds $\langle \sigma_{bs} \rangle$ for species $j$ at frequency
$i$, then solving the inverse problem means finding the vector of scatterer densities
$\mathbf{n} \ge 0$ which solves the equation 

$\mathbf{s}_v = \Sigma \mathbf{n}.$ (2)

The inverse approach has several advantages over the classify-then-integrate approach:
* It uses all available frequency information to integrate, in theory giving a more robust estimate of animal density than single-frequency integration.
* It extends naturally from a small number of narrowband frequencies to broadband spectra, or a mixture of the two.
* It can handle mixtures of different scatterers, a situation where dB differencing struggles.

However, even though the inverse approach has been around for some time, it has never
acheived widespread use in practice, because Equation 2 is often (if not usually)
underdetermined, and it involves inherent uncertainties that are nonlinear and 
potentially difficult to quantify.

### Bayesian Inverse Modeling

One way to address these challenges is to implement the inverse problem as a Bayesian
statistical model. Like all inversion methods, a Bayesian approach handles species mixtures, uses all 
available frequencies, and extends naturally to broadband signals. However, it has 
a few distinct advantages. The priors required for a Bayesian model provide a rigorous
way to incorporate assumptions, ecological knowledge, and/or information from direct 
sampling. Good priors provide the inverse model with additional information, improving
the quality of the solution and allowing even underdetermined  problems to be solved.
Bayesian models can incorporate multiple sources of uncertainty and propagate them
through to the solution, increasing the robustness ofthe results: a well-specified 
model should not produce solutions that are simultaneously wrong and confident.
Finally, these models are based on physical scattering processes, so their output is
fully interpretable, unlike some machine learning methods. Taken together, these 
advantages make the inverse approach robust and reliable enough to be used on real-world
acoustic data.

## Tutorial

### Installation

ProbabilisticEchoInversion.jl is a Julia package, so if you have not installed the Julia 
programming language, that's the first thing to do. You can download the latest version
it for free from the [official website](https://julialang.org/downloads/). Even better,
use the [Juliaup installer/version manager](https://github.com/JuliaLang/juliaup), which 
makes it much easier to upgrade when new Julia versions are released, and to maintain 
multiple Julia installations on your computer.

You can work with Julia files in any text editor, but for a nice integrated experience,
we can recommend [Visual Studio Code](https://code.visualstudio.com/) with the 
[Julia extension](https://www.julia-vscode.org/), [Jupyter](https://jupyter.org/), or 
[Pluto](https://plutojl.org/).

Once you have Julia installed, open the Julia command line (a.k.a. the read-evaluate-print-loop,
or REPL). While not required, it is easy and highly recommended to set up a local environment[^2] 
for each of your projects. To do that for this tutorial, run the following commands:

```julia
# create a new folder--could call it whatever you want
julia> mkdir("APESTutorial")

# change directory to the folder we just created
julia> cd("APESTutorial")

# type `]` to enter package manager mode
julia> ]

# activate the current directory as project environment
(@v1.8) pkg> activate .

(APESTutorial) pkg>
```

Install the package to this new project environment from GitHub by running the following
command:

```julia
(APESTutorial) pkg> add https://github.com/ElOceanografo/ProbabilisticEchoInversion.jl.git
```

Once it has downloaded and precompiled, you can exit the package manager by hitting backspace.
Load the package, and you're ready to go!

```julia
julia> using ProbabilisticEchoInversion
```
[^2] You don't strictly need to create a local environment, and can install 
ProbabilisticEchoInversion into the top-level Julia environment (i.e., `(@v1.8)` instead 
of `APESTutorial`). This will make it available automatically for all projects. However,
the more packages you install in the top-level environment, the more likely you are to end
up with conflicting versions and dependencies. In our experience, working with local
environments is *much* easier in the long run--and as a pleasant side effect, it makes it 
much easier to reproduce your analyses, since all the precise package versions you used
are recorded automatically in the Project.toml and Manifest.toml files.

### Loading and arranging your data

In a nutshell, you arrange your acoustic data in a multidimensional `DimArray` from 
DimensionalData.jl.[^1] For a typical downard-looking echosounder on a moving ship, your 
data would have three dimensions, depth x distance x frequency. You then specify the inverse 
problem you want to solve using the modeling language Turing.jl, and then use 
ProbabilisticEchoInversion to apply it to each depth/distance cell in your data.


[^1] (`DimArray`s behave like normal Julia arrays, but also 
let you index with named dimensions and have a bunch of nice functionality for subsetting
and plotting.) 

### Defining the model

### 

## More Advanced Examples

Because the inverse model is defined in Turing.jl's modeling language, APES is incredibly
flexible in terms of the data and situations to which it can be applied. For a more in-
depth look at some of its capabilities, please check out the fully-worked example 
problems in the [APESExamples repository](https://github.com/ElOceanografo/APESExamples).
These reproduce the analyses and figures in Urmy et al. 20**, and include:
* A simulated fish/zooplankton mixture, demonstrating how to solve simultaneously for scatterer size and density,
* A simulated mesopelagic scattering layer, demonstrating how to use ground-truth data to constrain the solution of an underdetermined inverse problem,
* An application of APES to a variety of mixed scattering types at the Aleutian shelfbreak
in the Gulf of Alaska, 
* And an application of APES to broadband backscatter from zooplankton and fish of mixed 
sizes in Barnabas Trough, south of Kodiak Island.


## Using, Citing, and Contributing
This software was developed by Sam Urmy at NOAA's Alaska Fisheries Science Center. As a 
product of the U.S. Government, it is free for anyone to use under a Creative Commons CC0
license. 

ProbabilisticEchoInversion.jl has been tested and peer-reviewed, but should still be
considered research-grade beta software rather than fully production-ready. If you try 
it on your own data, please do submit bug reports, comments, and feature requests via the
project's [GitHub repository](https://github.com/ElOceanografo/ProbabilisticEchoInversion.jl).
Pull requests are welcome, both for code and documentation!

Finally, if you do use APES in your own work, please cite the following publication:

TK

## API

```@index
```

```@autodocs
Modules = [ProbabilisticEchoInversion]
```
