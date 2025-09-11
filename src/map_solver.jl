"""
    MAPSolver([;optimizer, options])

Construct a `MAPSolver`, specifying how to invert a probabilistic backsattering
model using maximum a-posteriori optimization.  Default optimizer is LBFGS. See
Turing.jl and Optim.jl documentation for more information on available solvers
and options.
"""
Base.@kwdef struct MAPSolver <: AbstractSolver
    optimizer = LBFGS()
    args = ()
    kwargs = (;)
    hess_ad = AutoForwardDiff()
    verbose = false
end

struct MAPSolution{TM,TV,TO}
    mean::TM
    cov::TV
    optimizer::TO
end

function Base.show(io::IO, s::MAPSolution)
    print(io, s.optimizer.values)
end

function Base.show(io::IO, m::MIME"text/plain", s::MAPSolution)
    print(io, "MAPSolution with log-probability $(round(s.optimizer.lp, digits=2))")
    print(io, " and modal values:\n")
    show(io, m, s.optimizer.values)
end

Statistics.mean(s::MAPSolution) = s.mean
Statistics.cov(s::MAPSolution) = s.cov
function Statistics.cor(s::MAPSolution)
    D = diagm(1 ./ sqrt.(diag(s.cov)))
    return D * cov(s) * D
end
function Statistics.var(s::MAPSolution)
    v = diag(cov(s))
    v[v .< 0] .= NaN
    return v
end
Statistics.std(s::MAPSolution) = sqrt.(var(s))
cv(s; kwargs...) = std(s; kwargs...) ./ abs.(mean(s; kwargs...))
StatsBase.coef(s::MAPSolution) = coef(s.optimizer)

function solve(data, model::Function, solver::MAPSolver, params=())
    m = model(data, params)
    if solver.verbose 
        opt = maximum_a_posteriori(m, solver.optimizer, solver.args...; solver.kwargs...)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        opt = Logging.with_logger(logger) do
            maximum_a_posteriori(m, solver.optimizer, solver.args...; solver.kwargs...)
        end
        flush(io)
        close(io)
    end
    H = Symmetric(hessian(opt.f, solver.hess_ad, opt.optim_result.u))
    μ = opt.values
    μ_names = only(names(μ))
    C = try
        C = inv(H)
        NamedArray(C, names=(μ_names, μ_names))
    catch
        C = diagm(fill(Inf, size(H, 1)))
        NamedArray(C, names=(μ_names, μ_names))
    end
    return MAPSolution(μ, C, opt)
end