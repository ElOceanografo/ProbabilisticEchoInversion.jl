using Test
using ProbabilisticEchoInversion
using DimensionalData
using Statistics
using Turing

echogram = rand(X(0:0.1:2), Z(-10:1:0), F(30:42))

function mappointsslices(f, a, dims)
    dd = otherdims(a, dims)
    result = map(tup -> f(tup[1], a[tup[2]...]), zip(DimPoints(dd), DimKeys(dd)))
    return DimArray(result, dd)
end

f(coords, sv) = mean(sv)

post = mappointsslices(f, echogram, :F)


TS = randn(size(echogram, :F), 4)
prob = StaticAPESProblem(echogram[1, 1, :], TS, fill(Normal(), 4), Exponential(0.1))
sol = solve(prob, MCMCSolver())
