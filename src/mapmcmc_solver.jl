
"""
    MAPMCMCSolver([;optimizer, options])

Construct an `MAPMCMCSolver`, specifying how to invert a probabilistic backscattering
model using a combination of maximum a posteriori optimization and Markov-chain Monte 
Carlo. This simply means that an optimization routine finds the MAP point estimate of
the parameters, which is then used as the starting point for the MCMC run.

Arguments correspond exactly to the ones for `MAPSolver` and `MCMCSolver`; refer to 
their documentation for details.
"""
Base.@kwdef struct MAPMCMCSolver <: AbstractSolver
    sampler = NUTS(0.8)
    parallel = MCMCSerial()
    nsamples = 1000
    nchains = 1
    optimizer = LBFGS()
    args = ()
    kwargs = (progress=false,)
    verbose=false
end



function solve(data, model::Function, solver::MAPMCMCSolver, params=())
    m = model(data, params)
    if solver.verbose
        opt = maximum_a_posteriori(m, solver.optimizer, solver.args...; solver.kwargs...)
        return sample(m, solver.sampler, solver.parallel, solver.nsamples, solver.nchains;
            init_theta = opt.values.array, solver.kwargs...)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        chain = Logging.with_logger(logger) do
            opt = maximum_a_posteriori(m, solver.optimizer, solver.args...; solver.kwargs...)
            sample(m, solver.sampler, solver.parallel, solver.nsamples, solver.nchains;
                init_theta=opt.values.array, solver.kwargs...)
        end
        flush(io)
        close(io)
        return chain
    end
end