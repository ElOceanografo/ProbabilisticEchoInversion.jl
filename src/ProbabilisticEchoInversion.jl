module ProbabilisticEchoInversion

using Reexport
using CSV
using DataFrames, DataFramesMeta
using DimensionalData, DimensionalData.Dimensions
@reexport using Turing
using DifferentiationInterface
using OptimizationOptimJL
using NamedArrays
using ForwardDiff
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
    F,
    summarize_posterior

@dim F YDim "Frequency (kHz)"

abstract type AbstractSolver end
include("map_solver.jl")
include("mcmc_solver.jl")
include("mapmcmc_solver.jl")

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

_get_varnames(sol::Chains) = sol.name_map.parameters
_get_varnames(sol::MAPSolution) = only(names(sol.optimizer.values))
_get_varnames(sol::AbstractArray) = mapreduce(_get_varnames, union, sol)

function _get_group_varnames(sol, group)
    all_variables = _get_varnames(sol)
    return filter(v -> startswith(string(v), string(group)*"["), all_variables)
end

function _get_group_varnames(sol, groups::AbstractArray)
    group_variables = map(groups) do g
        _get_group_varnames(sol, g)
    end
    return vcat(group_variables...)
end

Base.getindex(sol::MAPSolution, key) = sol.optimizer.values[Symbol(key)]

_summarize(s::Chains, varname, f) = f(s[Symbol(varname)])
_summarize(s::MAPSolution, varname, f) = f(s)[Symbol(varname)]

regularize(sym) = Symbol(replace(string(sym),  "[" => "_", "]" => ""))

function summarize_posterior(f::Function, sol; variables=[], groups=[], regularize_names=true)
    variables = if isempty(variables) & isempty(groups)
        _get_varnames(sol) # default: all variables
    else
        variables
    end
    group_variables = _get_group_varnames(sol, groups)
    variables = union(variables, group_variables)
    summaries = map(variables) do varname
        f1 = s -> _summarize(s, varname, f)
        map(passmissing(f1), sol)
    end
    stack_names = regularize_names ? regularize.(variables) : variables
    return DimStack(summaries, name=stack_names)
end


end # module