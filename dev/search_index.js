var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ProbabilisticEchoInversion","category":"page"},{"location":"#ProbabilisticEchoInversion","page":"Home","title":"ProbabilisticEchoInversion","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ProbabilisticEchoInversion.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for ProbabilisticEchoInversion.jl!","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package is designed to solve the \"inverse problem\" for acoustic backscatter at multiple frequencies in a Bayesian framework. In other words, given observed echoes from one or more types of scatterers, this package will help you infer 1) what they were, ","category":"page"},{"location":"","page":"Home","title":"Home","text":"how many of them were present, and 3) how sure you can be about (1) and (2). We","category":"page"},{"location":"","page":"Home","title":"Home","text":"call this approach Automatic Probabilistic Echo Solving, or APES.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In a nutshell, you arrange your acoustic data in a multidimensional DimArray from  DimensionalData.jl.[1] For a typical downard-looking echosounder on a moving ship, your  data would have three dimensions, depth x distance x frequency.","category":"page"},{"location":"","page":"Home","title":"Home","text":"You then specify the inverse problem you want to solve using the modeling language Turing.jl,  and then use ProbabilisticEchoInversion to apply it to each depth/distance cell in your data","category":"page"},{"location":"","page":"Home","title":"Home","text":"[1] (DimArrays behave like normal Julia arrays, but also  let you index with named dimensions and have a bunch of nice functionality for subsetting and plotting.) ","category":"page"},{"location":"#Tutorial","page":"Home","title":"Tutorial","text":"","category":"section"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ProbabilisticEchoInversion]","category":"page"},{"location":"#ProbabilisticEchoInversion.MAPMCMCSolver","page":"Home","title":"ProbabilisticEchoInversion.MAPMCMCSolver","text":"MAPMCMCSolver([;optimizer, options])\n\nConstruct an MAPMCMCSolver, specifying how to invert a probabilistic backscattering model using a combination of maximum a posteriori optimization and Markov-chain Monte  Carlo. This simply means that an optimization routine finds the MAP point estimate of the parameters, which is then used as the starting point for the MCMC run.\n\nArguments correspond exactly to the ones for MAPSolver and MCMCSolver; refer to  their documentation for details.\n\n\n\n\n\n","category":"type"},{"location":"#ProbabilisticEchoInversion.MAPSolver","page":"Home","title":"ProbabilisticEchoInversion.MAPSolver","text":"MAPSolver([;optimizer, options])\n\nConstruct a MAPSolver, specifying how to invert a probabilistic backsattering model using maximum a-posteriori optimization.  Default optimizer is L-BFGS. See Turing.jl and Optim.jl documentation for more information on available solvers and options.\n\n\n\n\n\n","category":"type"},{"location":"#ProbabilisticEchoInversion.MCMCSolver","page":"Home","title":"ProbabilisticEchoInversion.MCMCSolver","text":"MCMCSolver([;sampler, parallel, nsamples, nchains; kwargs, verbose])\n\nConstruct an MCMCSolver, specifying how to invert a probabilistic backscattering model using Markov-chain Monte Carlo. By default uses the no-U-turn sampler with  acceptance rate 0.8 and collects 1000 samples. See Turing.jl documentation for more information on options for MCMC sampling.\n\n\n\n\n\n","category":"type"},{"location":"#ProbabilisticEchoInversion.apes-Tuple{DimensionalData.DimArray, Function, AbstractSolver}","page":"Home","title":"ProbabilisticEchoInversion.apes","text":"apes(echogram, model, solver[; params, result_handler, safe_precision, distributed])\n\nRun the automatic probabilistic echo solver defined by the inverse model and  solution method solver on the acoustic backstter data in echogram.\n\nArguments\n\nechogram::DimArray: Acoustic backscatter data in the linear domain (i.e.,   volume backsattering coefficient s_v, area backsattering coefficient s_a,   or nautical area scattering coefficient, NASC). The last dimension of    the DimArray should be named :F and index the acoustic frequencies;    all other dimensions should reference spatial/temporal coordinates.\nmodel::Function: Probabilistic inverse model defined with Turing.jl   or DynamicPPL.jl. This model should have the signature model(data, params),   where data and params contain the acoustic data and any additional   parameters. See below for more details.\nsolver::AbstractSolver: The method used to solve the inverse problem   specified in model. See MCMCSolver and MAPSolver for more detail.\nparams: Optional additional params to pass to model.\nresult_handler: Optional function to transform the output of the solver    before (for instance, by calculating the means of a Markov chain).\ndistributed::Bool=false: Whether to use all available processors when    fitting model to echogram cells.\n\nDetails\n\nThis function applies a probabilistic inverse backscattering model defined using Turing.jl to each spectrum in a mulifrequency echogram. The model's constructor must accept two arguments:\n\ndata: A NamedTuple or other structure accessible by dot-notation, with   fields coords, freqs, and backscatter. These contain the observed   acoustic data. The model doesn't have to use all of them.\nparams: Optional NamedTuple or other object, containing any constants or   auxiliary information used by the model.\n\n\n\n\n\n","category":"method"},{"location":"#ProbabilisticEchoInversion.iterspectra","page":"Home","title":"ProbabilisticEchoInversion.iterspectra","text":"iterspectra(echogram[, freqdim])\n\nGiven an mulifrequency or broadband echogram in the form of a DimArray, with the acoustic frequencies in one dimension (by default named :F), iterate over each spectrum. The iterator yields NamedTuples with three fields: \n\ncoords: Coordinates of the spectrum on the non-:F dimensions of the Array   (i.e. its location space/time)\nfreqs: Array of acoustic frequencies\nbackscatter: Array of backscatter values\n\n\n\n\n\n","category":"function"},{"location":"#ProbabilisticEchoInversion.mapspectra-Tuple{Any, DimensionalData.DimArray}","page":"Home","title":"ProbabilisticEchoInversion.mapspectra","text":"mapspectra(f, echogram[; freqdim, distributed])\n\nMap the function f over each spectrum in the DimArray echogram. By default assumes the acoustic frequencies are recorded in dimension :F, if this is not the case, specify the name of the dimension using the freqdim argument.\n\n\n\n\n\n","category":"method"},{"location":"#ProbabilisticEchoInversion.solve","page":"Home","title":"ProbabilisticEchoInversion.solve","text":"solve(data, model, solver[, params])\n\nRun the probabilistic inverse model defined by model on the acoustic backscatter spectrum in data, using solver as the inference engine.\n\n\n\n\n\n","category":"function"}]
}
