module ProbabilisticEchoInversion

using DimensionalData, DimensionalData.Dimensions
using Turing
using Optim
using Statistics, StatsBase
using Distributed
using DoubleFloats
using ProgressMeter
import Logging

export F,
    AbstractSolver,
    MCMCSolver,
    MAPSolver,
    solve,
    iterspectra,
    mapspectra,
    check_precision,
    apes

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
    kwargs = ()
    verbose = false
end


"""
    MAPSolver([;optimizer, options])

Construct a `MAPSolver`, specifying how to invert a probabilistic backsattering
model using maximum a-posteriori optimization.  Default optimizer is L-BFGS. See
Turing.jl and Optim.jl documentation for more information on available solvers
and options.
"""
Base.@kwdef struct MAPSolver <: AbstractSolver
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
        return sample(m, solver.sampler, solver.parallel, solver.nsamples, 
            solver.nchains; solver.kwargs...)
    end
    flush(io)
    close(io)
    return chain
end

function solve(data, model::Function, solver::MAPSolver, params=())
    m = model(data, params)
    if solver.verbose 
        return optimize(m, MAP(), solver.optimizer, solver.options)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        opt = optimize(m, MAP(), solver.optimizer, solver.options)
        flush(io)
        close(io)
    return opt
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
    pointskeys = zip(DimPoints(dd), DimKeys(dd))
    itr = ((coords=tup[1], freqs=freqs, backscatter=@view echogram[tup[2]...]) 
        for tup in pointskeys)
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

function check_precision(data, T=Double64)
    backscatter = data.backscatter
    if any(abs.(data.backscatter) .< sqrt(eps(eltype(data.backscatter))))
        bacscatter = T.(data.backscatter)
    end
    return (coords = data.coords, freqs=data.freqs, backscatter)
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
- `safe_precision::Bool=true`: Check whether backscatter absolute values are small
    enough to introduce floating-point errors. If so, convert them to higher-
    precision representation (default is `Double64`). Setting this to false can 
    speed up the calculations, but may introduce inference errors.
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
        params=(), result_handler=(x)->x, safe_precision=true, distributed=false)
    # numerical precision check here?
    # parallel options?
    function f(x)
        res = solve(x, model, solver, params)
        return result_handler(res)
    end
    return mapspectra(f, echogram, distributed=distributed)
end

end # module