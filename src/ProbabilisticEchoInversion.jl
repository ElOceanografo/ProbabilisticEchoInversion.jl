module ProbabilisticEchoInversion

using Reexport
using CSV
using DataFrames, DataFramesMeta
using DimensionalData, DimensionalData.Dimensions
@reexport using Turing
using Optim
using ForwardDiff
using FiniteDiff
using LinearAlgebra
using Statistics, StatsBase
using Distributed
using ProgressMeter
import Logging

export AbstractSolver,
    MCMCSolver,
    MAPSolver,
    MAPSolution,
    MAPMCMCSolver,
    solve,
    iterspectra,
    mapspectra,
    apes,
    cv,
    unstack_echogram,
    F

@dim F YDim "Frequency (kHz)"

abstract type AbstractSolver end

"""
    MCMCSolver([;sampler, parallel, nsamples, nchains; kwargs, verbose])

Construct an `MCMCSolver`, specifying how to invert a probabilistic backscattering
model using Markov-chain Monte Carlo. By default uses the no-U-turn sampler with 
acceptance rate 0.8 and collects 1000 samples. See Turing.jl documentation for more
information on options for MCMC sampling.
"""
Base.@kwdef struct MCMCSolver <: AbstractSolver
    sampler = NUTS(0.8)
    parallel = MCMCSerial()
    nsamples = 1000
    nchains = 1
    kwargs = (progress=false,)
    verbose = false
end


"""
    MAPSolver([;optimizer, options])

Construct a `MAPSolver`, specifying how to invert a probabilistic backsattering
model using maximum a-posteriori optimization.  Default optimizer is LBFGS. See
Turing.jl and Optim.jl documentation for more information on available solvers
and options.
"""
Base.@kwdef struct MAPSolver <: AbstractSolver
    optimizer = LBFGS()
    options = Optim.Options()
    hessian = :forwarddiff
    verbose=false

    function MAPSolver(optimizer, options, hessian, verbose)
        if hessian in [:forwarddiff, :finitediff]
            return new(optimizer, options, hessian, verbose)
        else
            throw(ArgumentError("`hessian` must be either :forwarddiff or :finitediff"))
        end
    end
end

"""
    MAPMCMCSolver([;optimizer, options])

Construct an `MAPMCMCSolver`, specifying how to invert a probabilistic backscattering
model using a combination of maximum a posteriori optimization and Markov-chain Monte 
Carlo. This simply means that an optimization routine finds the MAP point estimate of
the parameters, which is then used as the starting point for the MCMC run.

Arguments correspond exactly to the ones for `MAPSolver` and `MCMCSolver`; refer to 
their documentation for details.
"""
Base.@kwdef struct MAPMCMCSolver <: AbstractSolver
    sampler = NUTS(0.8)
    parallel = MCMCSerial()
    nsamples = 1000
    nchains = 1
    kwargs = (progress=false,)
    optimizer = LBFGS()
    options = Optim.Options()
    verbose=false
end


"""
    solve(data, model, solver[, params])

Run the probabilistic inverse model defined by `model` on the acoustic backscatter
spectrum in `data`, using `solver` as the inference engine.
"""
function solve(data, model::Function, solver::MCMCSolver, params=())
    m = model(data, params)
    if solver.verbose
        return sample(m, solver.sampler, solver.parallel, solver.nsamples, 
            solver.nchains; solver.kwargs...)
    end
    io = IOBuffer()
    logger = Logging.SimpleLogger(io, Logging.Error)
    chain = Logging.with_logger(logger) do
        sample(m, solver.sampler, solver.parallel, solver.nsamples, 
            solver.nchains; solver.kwargs...)
    end
    flush(io)
    close(io)
    return chain
end

function calculate_hessian(opt, solver)
    if solver.hessian == :forwarddiff
        return ForwardDiff.hessian(opt.f, opt.optim_result.minimizer)
    elseif solver.hessian == :finitediff
        return ForwardDiff.hessian(opt.f, opt.optim_result.minimizer)
    end
end

struct MAPSolution{TM,TV,TO}
    mean::TM
    cov::TV
    optimizer::TO
end

function Base.show(io::IO, s::MAPSolution)
    print(io, s.optimizer.values)
end

function Base.show(io::IO, m::MIME"text/plain", s::MAPSolution)
    print(io, "MAPSolution with log-probability $(round(s.optimizer.lp, digits=2))")
    print(io, " and modal values:\n")
    show(io, m, s.optimizer.values)
end

Statistics.mean(s::MAPSolution) = s.mean
Statistics.cov(s::MAPSolution) = s.cov
function Statistics.cor(s::MAPSolution)
    D = diagm(1 ./ sqrt.(diag(s.cov)))
    return D * cov(s) * D
end
Statistics.var(s::MAPSolution) = diag(cov(s))
Statistics.std(s::MAPSolution) = sqrt.(var(s))
cv(s; kwargs...) = std(s; kwargs...) ./ abs.(mean(s; kwargs...))
StatsBase.coef(s::MAPSolution) = coef(s.optimizer)

function solve(data, model::Function, solver::MAPSolver, params=())
    m = model(data, params)
    if solver.verbose 
        opt = optimize(m, MAP(), solver.optimizer, solver.options)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        opt = Logging.with_logger(logger) do
            optimize(m, MAP(), solver.optimizer, solver.options)
        end
        flush(io)
        close(io)
    end
    H = Symmetric(calculate_hessian(opt, solver))
    μ = opt.values
    try
        C = isposdef(H) ? inv(H) : pinv(H)
        return MAPSolution(μ, C, opt)
    catch
        C = diagm(fill(Inf, size(H, 1)))
        return MAPSolution(μ, C, opt)
    end
end

function solve(data, model::Function, solver::MAPMCMCSolver, params=())
    m = model(data, params)
    if solver.verbose
        opt = optimize(m, MAP(), solver.optimizer, solver.options)
        return sample(m, solver.sampler, solver.parallel, solver.nsamples, solver.nchains;
            init_theta = opt.values.array, solver.kwargs...)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        chain = Logging.with_logger(logger) do
            opt = optimize(m, MAP(), solver.optimizer, solver.options)
            sample(m, solver.sampler, solver.parallel, solver.nsamples, solver.nchains;
                init_theta=opt.values.array, solver.kwargs...)
        end
        flush(io)
        close(io)
        return chain
    end
end

"""
    iterspectra(echogram[, freqdim])

Given an mulifrequency or broadband echogram in the form of a `DimArray`, with the
acoustic frequencies in one dimension (by default named `:F`), iterate over each
spectrum. The iterator yields `NamedTuples` with three fields: 
- `coords`: Coordinates of the spectrum on the non-`:F` dimensions of the Array
    (i.e. its location space/time)
- `freqs`: Array of acoustic frequencies
- `backscatter`: Array of backscatter values
"""
function iterspectra(echogram::DimArray, freqdim=:F)
    freqs = collect(dims(echogram, freqdim))
    dd = otherdims(echogram, freqdim)
    ddnames = name.(dd)
    pointskeys = zip(DimPoints(dd), DimKeys(dd))
    itr = (
        (
            coords=NamedTuple(zip(ddnames, tup[1])), 
            freqs=freqs,
            backscatter=vec(@view echogram[tup[2]...])
        )
        for tup in pointskeys
    )
    return itr
end

"""
    mapspectra(f, echogram[; freqdim, distributed])

Map the function `f` over each spectrum in the `DimArray` `echogram`. By default
assumes the acoustic frequencies are recorded in dimension `:F`, if this is not
the case, specify the name of the dimension using the `freqdim` argument.
"""
function mapspectra(f, echogram::DimArray; freqdim=:F, distributed=false)
    itr = iterspectra(echogram, freqdim)
    if distributed
        result = @showprogress pmap(f, itr)
    else
        result = @showprogress map(f, itr)
    end
    dd = otherdims(echogram, freqdim)
    return DimArray(result, dd)
end

"""
    apes(echogram, model, solver[; params, result_handler, safe_precision, distributed])

Run the automatic probabilistic echo solver defined by the inverse `model` and 
solution method `solver` on the acoustic backstter data in `echogram`.

# Arguments
- `echogram::DimArray`: Acoustic backscatter data in the linear domain (i.e.,
    volume backsattering coefficient ``s_v``, area backsattering coefficient ``s_a``,
    or nautical area scattering coefficient, ``NASC``). The last dimension of 
    the `DimArray` should be named `:F` and index the acoustic frequencies; 
    all other dimensions should reference spatial/temporal coordinates.
- `model::Function`: Probabilistic inverse model defined with Turing.jl
    or DynamicPPL.jl. This model should have the signature `model(data, params)`,
    where `data` and `params` contain the acoustic data and any additional
    parameters. See below for more details.
- `solver::AbstractSolver`: The method used to solve the inverse problem
    specified in `model`. See `MCMCSolver` and `MAPSolver` for more detail.
- `params`: Optional additional params to pass to `model`.
- `result_handler`: Optional function to transform the output of the solver 
    before (for instance, by calculating the means of a Markov chain).
- `distributed::Bool=false`: Whether to use all available processors when 
    fitting model to echogram cells.

# Details

This function applies a probabilistic inverse backscattering model defined using
Turing.jl to each spectrum in a mulifrequency echogram. The model's constructor
must accept two arguments:
- `data`: A `NamedTuple` or other structure accessible by dot-notation, with
    fields `coords`, `freqs`, and `backscatter`. These contain the observed
    acoustic data. The model doesn't have to use all of them.
- `params`: Optional `NamedTuple` or other object, containing any constants or
    auxiliary information used by the model.
"""
function apes(echogram::DimArray, model::Function, solver::AbstractSolver; 
        params=(), result_handler=(x)->x, distributed=false)
    function f(x)
        if all(ismissing.(x.backscatter))
            return missing
        else
            res = solve(x, model, solver, params)
            return result_handler(res)
        end
    end
    return mapspectra(f, echogram, distributed=distributed)
end

"""
    unstack_echogram(echo_df, xcol, ycol, fcol, svcol[, X, Y, F])

Convert a long-format `DataFrame` of multifrequency acoustic data into a three-dimensional
`DimArray` representing an echogram. This `DataFrame` must contain the following columns:

1. An x-coordinate, such as along-track distance or time
2. A y-coordinate, such as depth or range from the transducer
3. Acoustic frequency, and
4. Acoustic mean volume backscatter (can be in linear or decibel units )

The names of these columns are passed to `unstack_echogram` as the second through
fifth arguments, respectively.

- `echo_df::DataFrame` : Long-format `DataFrame` with columns for an x and y coordinate,
acoustic frequency, and backscatter.
- `xcol`, `ycol`, `fcol`, `svcol` : Names of the columns in `echo_df` (as `Symbol`s)
- `X`, `Y`, `F` : Optional `Dimensions` to apply to the x, y, and frequency axes of the 
resulting `DimArray`.

"""
function unstack_echogram(echo_df::DataFrame, xcol::Symbol, ycol::Symbol, fcol::Symbol, svcol::Symbol, 
    # ::Type{X}, ::Type{Y}, ::Type{F}) where {X<:Dimension, Y<:Dimension, F<:Dimension}
    X=X, Y=Y, F=F)
    freqs = sort(unique(echo_df[:, fcol]))
    x = sort(unique(echo_df[:, xcol]))
    y = sort(unique(echo_df[:, ycol]))
    stack = map(freqs) do f 
        @chain echo_df begin
            subset(fcol => x -> x .== f)
            sort([xcol, ycol])
            unstack(ycol, xcol, svcol)
            sort(ycol)
            select(Not(ycol))
            Array()
        end
    end
    echogram = cat(stack..., dims=3)
    return DimArray(echogram, (Y(y), X(x), F(freqs)))
end
end # module