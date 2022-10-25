using Distributed
addprocs(exeflags="--project=$(Base.active_project())")

using ProbabilisticEchoInversion
@everywhere begin
    # using Pkg; Pkg.activate(".")
    using Test
    using ProbabilisticEchoInversion
    using DimensionalData
    using Statistics
    using Turing
end

echogram = rand(X(0:0.1:2), Z(-10:1:0), F(30:42))

iterspectra(echogram)

@everywhere f(data) = data.freqs' * data.backscatter / sum(data.backscatter)
map(f, iterspectra(echogram))
mapspectra(f, echogram)
mapspectra(f, echogram, distributed=true)

TS = randn(size(echogram, :F), 4)

@everywhere @model function examplemodel(data, params)
    nspp = size(params.TS, 2)
    Σ = exp10.(params.TS / 10)
    logn ~ filldist(Normal(0, 1), nspp)
    n = exp.(logn)
    η ~ Exponential(0.1)
    ϵ = data.backscatter .* η
    data.backscatter .~ Normal.(Σ * n, ϵ)
end

par = (TS = TS,)


data = first(iterspectra(echogram))
mod = examplemodel(data, par)
solver = MCMCSolver(nsamples=100, kwargs=(progress=false,))

solve(data, examplemodel, MCMCSolver(verbose=true), par)
solve(data, examplemodel, MCMCSolver(parallel=MCMCThreads(), nchains=4), par)
solve(data, examplemodel, solver, par)

solve(data, examplemodel, MAPSolver(), par)

result = apes(echogram, examplemodel, solver, params=par, result_handler=mean)
result = apes(echogram, examplemodel, solver, params=par, result_handler=mean, 
    safe_precision=false, distributed=true)
result = apes(echogram, examplemodel, MAPSolver(), params=par)
result = apes(echogram, examplemodel, MAPSolver(), params=par, distributed=true)

@test size(result, 1) == size(echogram, 1)
@test size(result, 2) == size(echogram, 2)