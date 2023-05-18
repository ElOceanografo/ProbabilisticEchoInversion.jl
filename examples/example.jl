using ProbabilisticEchoInversion
using CSV, DataFrames, DataFramesMeta
using DimensionalData, DimensionalData.Dimensions
using Plots

freqs = [18, 38, 70, 120, 200]

echo_df = map(freqs) do f
    filename = joinpath(@__DIR__, "DY1702_haul24_$(f)kHz.csv")
    df = CSV.read(filename, DataFrame)
    df[!, :frequency] .= f
    return df
end
echo_df = vcat(echo_df...)

@dim Z YDim "Depth (m)"
@dim D XDim "Distance (km)"
echo = unstack_echogram(echo_df, :Dist_M, :Layer_depth_min, :frequency, :Sv_mean, D, Z, F)

heatmap(echo[F(At(120))], yflip=true)
savefig(joinpath(@__DIR__, "echogram.png"))

@model function examplemodel(data, params)
    nfreq, nspp = size(params.TS)
    Σ = exp10.(params.TS ./ 10)
    # priors
    logn ~ arraydist(Normal.(zeros(nspp), fill(3, nspp)))
    ϵ ~ Exponential(1.0)
    # Predict Sv based on scatterer density and TS
    n = exp10.(logn)
    μ = 10log10.(Σ * n)
    # Compare observed to predicted
    data.backscatter .~ Normal.(μ, fill(ϵ, nfreq))
end


TS_shrimp = [-100, -90, -82, -76.2, -73.7]
TS_pollock = [-34.6, -35.0, -35.6, -36.6, -38.5]
TS_myctophid = [-73.0, -58.0, -65, -67.2, -70.0]
TS = [TS_pollock TS_myctophid TS_shrimp]
params = (; TS)

solution_mcmc = apes(echo, examplemodel, MCMCSolver(), params=params);
post_mean = passmissing(mean).(solution_mcmc);
post_cv = passmissing(cv).(solution_mcmc);

solution_map = apes(echo, examplemodel, MAPSolver(), params=params);
post_mean = passmissing(mean).(solution_map);
post_cv = passmissing(cv).(solution_map);


