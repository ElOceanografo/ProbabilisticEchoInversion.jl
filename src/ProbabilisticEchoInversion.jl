module ProbabilisticEchoInversion

using DimensionalData, DimensionalData.Dimensions
using Turing
using Statistics

export F,
    StaticAPESProblem,
    StaticAPESModel,
    MCMCSolver,
    solve

@dim F YDim "Frequency (kHz)"


# echomodel(sv, nprior, metadata, params)


@model function StaticAPESModel(sv, TS, logn_prior, cv_prior, sv_inflated, sv_inflation)
    nfreqs = size(TS, 1)
    Σbs = exp10.(TS ./ 10)
    cv ~ cv_prior
    ηsv = max.(cv .* sv, sqrt(eps()) * sv_inflation)
    logn ~ arraydist(logn_prior)
    n = exp.(logn)
    sv_pred = Σbs * n * sv_inflation

    for i in 1:nfreqs
        if ! ismissing(sv_inflated)
            sv_inflated[i] ~ Normal(sv_pred[i], ηsv[i])
        end
    end
end

struct StaticAPESProblem
    sv
    sv_inflation
    sv_inflated
    TS
    n_prior
    cv_prior
    model
end

function StaticAPESProblem(sv, TS, n_prior, cv_prior)
    sv_inflation = all(ismissing.(sv)) ? 1.0 : 1/maximum(skipmissing(sv))
    sv_inflated = sv .* sv_inflation
    model = StaticAPESModel(sv, TS, n_prior, cv_prior, sv_inflated, sv_inflation)
    return StaticAPESProblem(
        sv, 
        sv_inflation, 
        sv_inflated, 
        TS, 
        n_prior, 
        cv_prior, 
        model)
end

struct MCMCSolver
    sampler
    nsamples
    args
    kwargs
    summarizer
end
MCMCSolver() = MCMCSolver(NUTS(0.8), 1000, (), (), x -> x)

function solve(prob::StaticAPESProblem, solver::MCMCSolver)
    return sample(prob.model, solver.sampler, solver.nsamples, solver.args..., solver.kwargs...)
end

end # module