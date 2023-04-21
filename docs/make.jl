using ProbabilisticEchoInversion
using Documenter

DocMeta.setdocmeta!(ProbabilisticEchoInversion, :DocTestSetup, :(using ProbabilisticEchoInversion); recursive=true)

makedocs(;
    modules=[ProbabilisticEchoInversion],
    authors="Sam Urmy <oceanographerschoice@gmail.com> and contributors",
    repo="https://github.com/user/ProbabilisticEchoInversion.jl/blob/{commit}{path}#{line}",
    sitename="ProbabilisticEchoInversion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://user.github.io/ProbabilisticEchoInversion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/user/ProbabilisticEchoInversion.jl",
    devbranch="main",
)
