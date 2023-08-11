using ProbabilisticEchoInversion
using Documenter

DocMeta.setdocmeta!(ProbabilisticEchoInversion, :DocTestSetup, :(using ProbabilisticEchoInversion); recursive=true)

makedocs(;
    modules=[ProbabilisticEchoInversion],
    authors="Sam Urmy <oceanographerschoice@gmail.com> and contributors",
    repo="https://github.com/ElOceanografo/ProbabilisticEchoInversion.jl/blob/{commit}{path}#{line}",
    sitename="ProbabilisticEchoInversion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ElOceanografo.github.io/ProbabilisticEchoInversion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ElOceanografo/ProbabilisticEchoInversion.jl",
    devbranch="main",
    versions=["stable" => "v^"],
)
