using Distributed
addprocs(exeflags="--project=$(Base.active_project())")

using ProbabilisticEchoInversion
using DataFrames
using Random

@everywhere begin
    # using Pkg; Pkg.activate(".")
    using Test
    using ProbabilisticEchoInversion
    using DimensionalData
    using Statistics
    using Turing
end

Random.seed!(12345)

echo_df = allcombinations(DataFrame, X=0:0.1:2, Z=-10:1:0, F=30:42)
echo_df.backscatter .= rand.()
echogram = unstack_echogram(echo_df, :X, :Z, :F, :backscatter)

cell_match = map(eachrow(echo_df)) do row
    echogram[X(At(row.X)), Y(At(row.Z)), F(At(row.F))] == row.backscatter
end
@test all(cell_match)

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
solve(data, examplemodel, MAPMCMCSolver(verbose=true), par)

result = apes(echogram, examplemodel, solver, params=par, result_handler=mean)
result = apes(echogram, examplemodel, solver, params=par, result_handler=mean,
    distributed=true)
result = apes(echogram, examplemodel, MAPSolver(), params=par)
result = apes(echogram, examplemodel, MAPSolver(), params=par, distributed=true)

@test size(result, 1) == size(echogram, 1)
@test size(result, 2) == size(echogram, 2)

for r in result
    @test mean(r) == r.mean
    @test cov(r) == r.cov
    @test mean(r) == r.optimizer.values
    var(r)
    std(r)
    cv(r)
end

