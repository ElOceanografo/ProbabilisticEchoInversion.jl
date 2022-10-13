module ProbabilisticEchoInversion

using DimensionalData, DimensionalData.Dimensions
using Turing
using Statistics
using DoubleFloats
using ProgressMeter
import Logging

export F,
    AbstractSolver,
    MCMCSolver,
    solve,
    iterspectra,
    mapspectra,
    check_precision,
    apes

@dim F YDim "Frequency (kHz)"

abstract type AbstractSolver end

Base.@kwdef struct MCMCSolver <: AbstractSolver
    sampler = NUTS(0.8)
    nsamples = 1000
    args = ()
    kwargs = ()
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
        return sample(m, solver.sampler, solver.nsamples, solver.args...; solver.kwargs...)
    end
    logger = Logging.SimpleLogger(Logging.Error)
    chain = Logging.with_logger(logger) do
        sample(m, solver.sampler, solver.nsamples, solver.args...; solver.kwargs...)
    end
    return chain
end

function iterspectra(echogram, freqdim=:F)
    freqs = collect(dims(echogram, freqdim))
    dd = otherdims(echogram, freqdim)
    pointskeys = zip(DimPoints(dd), DimKeys(dd))
    itr = ((coords=tup[1], freqs=freqs, backscatter=@view echogram[tup[2]...]) for tup in pointskeys)
    return itr
end

function mapspectra(f, echogram, freqdim=:F)
    itr = iterspectra(echogram, freqdim)
    result = @showprogress map(f, itr)
    dd = otherdims(echogram, freqdim)
    return DimArray(result, dd)
end

function check_precision(data, t=Double64)
    if any(abs.(data.backscatter) .< sqrt(eps(eltype(data.backscatter))))
        bacscatter = t.(data.backscatter)
    else
        backscatter = data.backscatter
    end
    return (coords = data.coords, freqs=data.freqs, backscatter)
end


"""
    apes(echogram, model, solver[; params, result_handler, safe_precision])

Run the automatic probabilistic echo solver defined by the inverse `model` and 
solution method `solver` on the acoustic backstter data in `echogram`.

# Arguments
- `echogram::DimArray`: Acoustic backscatter data in the linear domain. The 
    last dimension of the `DimArray` should be named `:F` and index the
    acoustic frequencies; all other dimensions should reference spatial/temporal
    coordinates.
- `model::Function`: Probabilistic inverse model defined with Turing.jl
    or DynamicPPL.jl. This model should have the signature `model(data, params)`,
    where `data` and `params` are `NamedTuple`s containing acoustic data and any
    additional parameters. See below for more details.
- `solver::AbstractSolver`: The method used to solve the inverse problem
    specified in `model`. See `MCMCSolver` and `MAPSolver` for more detail.
- `params`: Optional additional params to pass to `model`.
- `result_handler`: Optional function to transform the output of the solver 
    before (for instance, by calculating the means of a Markov chain).
- `safe_precision::Bool=true`: Check whether backscatter absolute values are small
    enough to introduce floating-point errors. If so, convert them to higher-
    precision representation (default is `Double64`). Setting this to false may 
    speed up inference, but 

# Details


"""
function apes(echogram::DimArray, model::Function, solver::AbstractSolver; 
        params=(), result_handler=(x)->x, safe_precision=true)
    # numerical precision check here?
    # parallel options?
    function f(x)
        x = safe_precision ? check_precision(x) : x
        res = solve(x, model, solver, params)
        return result_handler(res)
    end
    return mapspectra(f, echogram)
end

end # module