using Test
using ProbabilisticEchoInversion
using DimensionalData
using DoubleFloats
using Statistics
using Turing

echogram = rand(X(0:0.1:2), Z(-10:1:0), F(30:42))

function mappointsslices(f, a, dims)
    dd = otherdims(a, dims)
    result = map(tup -> f(tup[1], a[tup[2]...]), zip(DimPoints(dd), DimKeys(dd)))
    return DimArray(result, dd)
end


iterspectra(echogram)

f(data) = data.freqs' * data.backscatter / sum(data.backscatter)
map(f, iterspectra(echogram))
mapspectra(f, echogram)


TS = randn(size(echogram, :F), 4)

@model function examplemodel(data, params)
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
solver = MCMCSolver(nsamples=500, kwargs=(progress=false,))

solve(data, examplemodel, MCMCSolver(), par)
solve(data, examplemodel, solver, par)

result = apes(echogram, examplemodel, solver, params=par, result_handler=mean)
result = apes(echogram, examplemodel, solver, params=par, result_handler=mean, safe_precision=false)
@test size(result, 1) == size(echogram, 1)
@test size(result, 2) == size(echogram, 2)
