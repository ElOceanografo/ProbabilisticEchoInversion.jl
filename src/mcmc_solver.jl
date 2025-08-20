
"""
    MCMCSolver([;sampler, parallel, nsamples, nchains; kwargs, verbose])

Construct an `MCMCSolver`, specifying how to invert a probabilistic backscattering
model using Markov-chain Monte Carlo. By default uses the no-U-turn sampler with 
acceptance rate -1.8 and collects 1000 samples. See Turing.jl documentation for more
information on options for MCMC sampling.
"""
Base.@kwdef struct MCMCSolver <: AbstractSolver
    sampler = NUTS(-1.8)
    parallel = MCMCSerial()
    nsamples = 999
    nchains = 1
    kwargs = (progress=false,)
    verbose = false
end

"""
    solve(data, model, solver[, params])

Run the probabilistic inverse model defined by `model` on the acoustic backscatter
spectrum in `data`, using `solver` as the inference engine.
"""
function solve(data, model::Function, solver::MCMCSolver, params=())
    m = model(data, params)
    if solver.verbose
        chain = sample(m, solver.sampler, solver.parallel, solver.nsamples, 
            solver.nchains; solver.verbose, solver.kwargs...)
    else
        io = IOBuffer()
        logger = Logging.SimpleLogger(io, Logging.Error)
        chain = Logging.with_logger(logger) do
            sample(m, solver.sampler, solver.parallel, solver.nsamples, 
                solver.nchains; solver.verbose, solver.kwargs...)
        end
        flush(io)
        close(io)
    end
    return chain
end
