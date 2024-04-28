using ProbabilisticEchoInversion
using CSV, DataFrames 
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
echo = unstack_echogram(echo_df, :Dist_M, :Layer_depth_min, :frequency, :Sv_mean, D, Z)

heatmap(echo[F(At(120))], yflip=true)
savefig(joinpath(@__DIR__, "echogram.png"))

@model function examplemodel(data, params)
    nfreq, nspp = size(params.TS)
    Σ = exp10.(params.TS ./ 10)

    # define priors
    logn ~ arraydist(Normal.(zeros(nspp), fill(3, nspp))) # scatterer log-densities
    ϵ ~ Exponential(1.0)

    # Predict Sv based on scatterer density and TS
    n = exp10.(logn)
    μ = 10log10.(Σ * n)

    # Compare observed to predicted backscatter
    data.backscatter .~ Normal.(μ, fill(ϵ, nfreq))
end


TS_shrimp = [-100, -90, -82, -76.2, -73.7]
TS_pollock = [-34.6, -35.0, -35.6, -36.6, -38.5]
TS_myctophid = [-73.0, -58.0, -65, -67.2, -70.0]
TS = [TS_pollock TS_myctophid TS_shrimp]
params = (; TS)

plot(freqs, TS, marker=:o, label=["Pollock" "Myctophid" "Shrimp"],
    xlabel="Frequency (kHz)", ylabel="TS (dB re m⁻²)")
savefig(joinpath(@__DIR__, "ts.png"))

solution_mcmc = apes(echo, examplemodel, MCMCSolver(), params=params);

function posterior_statistic(stat, solution, var)
    f = passmissing(chn -> mapslices(stat, Array(group(chn, var)), dims=1))
    return f.(solution)
end
pm = posterior_statistic(mean, solution_mcmc, :logn);
vec(pm)
size(pm)
pm = reduce(hcat, vec(pm))
replace(pm, missing)

post_mean = passmissing(chn -> mean(Array(chn), dims=1)).(solution_mcmc);
post_cv = passmissing(chn -> cv(Array(chn), dims=1)).(solution_mcmc);

post_pollock = map(passmissing(m -> m[1]), post_mean)
post_myctophid = map(passmissing(m -> m[2]), post_mean)
post_shrimp = map(passmissing(m -> m[3]), post_mean)

post_cv_pollock = map(passmissing(m -> m[1]), post_cv)
post_cv_myctophid = map(passmissing(m -> m[2]), post_cv)
post_cv_shrimp = map(passmissing(m -> m[3]), post_cv)

plot(
    heatmap(post_pollock),
    heatmap(post_myctophid),
    heatmap(post_shrimp),
    yflip=true, clims=(-6, 3)
)
plot(
    heatmap(post_cv_pollock, c=:viridis),
    heatmap(post_cv_myctophid, c=:viridis),
    heatmap(post_cv_shrimp, c=:viridis),
    yflip=true, clims=(0, 10)
)


solution_map = apes(echo, examplemodel, MAPSolver(), params=params);
post_mean_map = passmissing(mean).(solution_map);
post_cv_map = passmissing(cv).(solution_map);

