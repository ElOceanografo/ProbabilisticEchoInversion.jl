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

### Installation and setup

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
or REPL). While not required, it is easy and highly recommended to set up a local environment[^1] 
for each of your projects. To do that for this tutorial, run the following commands:

```julia
# create a new folder--could call it whatever you want
julia> mkdir("APESTutorial")

# change directory to the folder we just created
julia> cd("APESTutorial")

# type `]` to enter package manager mode
julia> ]

# activate the current directory as project environment
(@v1.9) pkg> activate .

(APESTutorial) pkg>
```

Install the package to this new project environment from GitHub by running the following
command:

```julia
(APESTutorial) pkg> add https://github.com/ElOceanografo/ProbabilisticEchoInversion.jl.git
```

Once it has downloaded and precompiled, you can exit the package manager by hitting backspace.
To run the rest of this tutorial yourself, you'll need the data files located
[here](https://github.com/ElOceanografo/ProbabilisticEchoInversion.jl/tree/main/examples).
Download them to the project directory you just created. You can also download the `example.jl`
script, which contains all the following code in one place.

[^1] You don't strictly need to create a local environment, and can install 
ProbabilisticEchoInversion into the top-level Julia environment (i.e., `(@v1.9)` instead 
of `APESTutorial`). This will make it available automatically for all projects. However,
the more packages you install in the top-level environment, the more likely you are to end
up with conflicting versions and dependencies. In our experience, working with local
environments is *much* easier in the long run--and as a pleasant side effect, it makes it 
much easier to reproduce your analyses, since all the precise package versions you used
are recorded automatically in the Project.toml and Manifest.toml files.

### Loading and arranging your data

ProbabilisticEchoInversion expects your acoustic data to be arranged in a multidimensional
`DimArray` from DimensionalData.jl. `DimArray`s behave like normal Julia arrays, but also 
let you index with named dimensions and have a bunch of nice functionality for subsetting,
slicing, and plotting. For a typical downard-looking echosounder on a moving ship, your 
data have three dimensions - depth, distance, and frequency - meaning this array will be 
three-dimensional.

If your data are already stored this way, for instance in a NetCDF or .mat file, it will
be easy to convert them to a `DimArray` (refer to the
[DimensionalData.jl documentation](https://rafaqz.github.io/DimensionalData.jl/stable/course/)
for how to do that). If your data are stored as a table in "long" format, as is typical
for .CSV exports from Echoview, you will need to do some reshaping first. This package
provides a function `unstack_echogram` to do that for you.

First, load the required packages.
```julia
using ProbabilisticEchoInversion
using CSV, DataFrames
using DimensionalData, DimensionalData.Dimensions
```

There are five comma-delimited data files, one for each frequency. We read them in, add a
`frequency` column to each one, and use `vcat` to stack them all into a single `DataFrame`.
```julia
freqs = [18, 38, 70, 120, 200]
echo_df = map(freqs) do f
    filename = joinpath(@__DIR__, "DY1702_haul24_$(f)kHz.csv")
    df = CSV.read(filename, DataFrame)
    df[!, :frequency] .= f
    return df
end
echo_df = vcat(echo_df...)
```
Next, we'll transform this data frame into a 3-d `DimArray`. The `unstack_echogram`
function takes long-format a `DataFrame` as its first argument. That `DataFrame`
needs to have at least four columns:

1. An x-coordinate, such as along-track distance or time
2. A y-coordinate, such as depth or range from the transducer
3. Acoustic frequency, and
4. Acoustic mean volume backscatter (can be in linear or decibel units )

The names of these columns are passed to `unstack_echogram` as the second through
fifth arguments, respectively. In this example, the data files are standard .csv
exports from Echoview, and the relevant column names are:
```julia
echo = unstack_echogram(echo_df, :Dist_M, :Layer_depth_min, :frequency, :Sv_mean)
```
By default, the axes of `echo` will be named `X`, `Y`, and `F`. If you'd like to 
define your own dimensions, this is easy to do with the `@dim` macro from
`DimensionalData.Dimensions`. You can then supply them in the optional final
three arguments to `unstack_echogram` and they will be applied to the 
`DimArray` it returns.
```julia
@dim Z YDim "Depth (m)"
@dim D XDim "Distance (km)"
echo = unstack_echogram(echo_df, :Dist_M, :Layer_depth_min, :frequency, :Sv_mean, D, Z)
```
It is now easy to manipulate the multifrequency echogram, for instance by selecting a 
slice by frequency and plotting it. Refer to the DimensionalData.jl docs to learn more
about how to slice and dice `DimArrays`.
```julia
heatmap(echo[F(At(120))], yflip=true)
```

### Defining the model

Once the data are loaded in, we need to define the inverse model we want to solve.
This is done using the probabilistic programming language Turing.jl. If you are familiar
with BUGS, JAGS, or Stan, model definitions in Turing are conceptually very similar. If 
you have not worked with it before, it is worth studying the Turing 
[documentation](https://turinglang.org/dev/docs/using-turing/get-started) before going
any further.

A very simple inverse model is defined below.

```julia
@model function examplemodel(data, params)

    nfreq, nspp = size(params.TS)
    Σ = exp10.(params.TS ./ 10)

    # define priors
    logn ~ arraydist(Normal.(zeros(nspp), fill(3, nspp))) # scatterer log-densities
    ϵ ~ Exponential(1.0) # observation error variance

    # Predict Sv based on scatterer density and TS
    n = exp10.(logn)
    μ = 10log10.(Σ * n)

    # Compare observed to predicted backscatter
    data.backscatter .~ Normal.(μ, fill(ϵ, nfreq))
end

```
To work with APES, your model function must accept two arguments. The first, `data`, 
contains the observed acoustic data. It will be a `NamedTuple` with fields named `coords`,
`freqs`, and `backscatter` that gets generated automatically for each acoustic cell When
you call the `apes` function to actually run the analysis. You can use any of these
fields inside the model if you want, but `data.backscatter` is the most important, since
it contains the actual observations.

The second argument, `params` containins any constants or auxiliary information you want to
pass to the model. It will typically be a `NamedTuple`, but can be any type of object. If
your model doesn't need any other info, you can just supply an empty tuple `()`.
Here, `params` is going to hold a single item, a matrix of target strengths (TS).

This model assumes a fixed number of scattering classes are present, each with a known
TS spectrum. It puts a vague prior on their log-densities, and assumes a single error
variance for all frequencies. 

> ⚠ Note that this model is defined in the logarithmic domain - that is, the scatterer
> densities are written as log-densities, and the observed data are assumed to be 
> decibel-valued mean volume backscattering strengths ($S_v$) instead of linear mean
> volume backscattering coefficients ($s_v$). While not strictly required, defining
> your models this way is a *really good idea*. The small absolute values and wide
> ranges of both scatterer densities and observed backscatter means that linear-domain
> models often have problems with floating-point precision that can manifest in
> inefficient and/or incorrect inference.

The last task to set up our model is to choose our candidate scatterers and set up
the TS matrix we are going to pass to the model via `params`. A research trawl performed
at this location found a mixture of Alaska pollock (*Gadus chalcogrammus*), unidentified
lanternfish, and Pacific glass shrimp (*Pasiphaea pacifica*). We will assume these were
the main scatterers present and define three TS spectra at our five frequencies. We then 
concatenate them into a matrix and pack it into our named tuple of parameters.

```julia
TS_shrimp = [-100, -90, -82, -76.2, -73.7]
TS_pollock = [-34.6, -35.0, -35.6, -36.6, -38.5]
TS_myctophid = [-73.0, -58.0, -65, -67.2, -70.0]
TS = [TS_pollock TS_myctophid TS_shrimp]
params = (; TS)
```

### Running the model

Once the data, parameters, and model are all set up, running it is just one line of code.

```julia
solution_mcmc = apes(echo, examplemodel, MCMCSolver(), params=params)
```
This will draw 1,000 samples from the joint posterior of the model for each acoustic cell,
using the No-U-Turn Sampler (NUTS) from Turing. Any cells where all backscatter values
are `missing` (e.g., below-bottom cells) will be skipped. Altogether, the inference will
take on the order of 5-20 minutes to finish, depending on your machine.

If you don't have time to wait for (asymptotically) exact inference, you can opt for a much
faster option: maximum-a-posteriori optimization, with errors estimated via the delta method.
This is done by changing the solver argument to `MAPSolver()`. Inference in this case
takes just a few seconds.

```julia
solution_map = apes(echo, examplemodel, MapSolver(), params=params)
```

In either case (MCMC chains or optimization fits) the results are returned in a `DimArray`
that shares the first two dimensions as `echo`, so they are also easy to manipulate and 
plot. For instance, arrays of posterior means and standard deviations can be obtained this way,

```julia
post_mean = passmissing(mean).(solution_mcmc)
post_mean = passmissing(std).(solution_mcmc)
```

where we use `passmissing` to deal with the fact that some of the result cells contain
MCMC chains and some are missing.

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
