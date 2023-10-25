using PMSimulatorCallbacks
using Documenter

DocMeta.setdocmeta!(PMSimulatorCallbacks, :DocTestSetup, :(using PMSimulatorCallbacks); recursive=true)

makedocs(;
    modules=[PMSimulatorCallbacks],
    authors="Timothy Knab",
    repo="https://github.com/timknab/PMSimulatorCallbacks.jl/blob/{commit}{path}#{line}",
    sitename="PMSimulatorCallbacks.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://timknab.github.io/PMSimulatorCallbacks.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/timknab/PMSimulatorCallbacks.jl",
    devbranch="main",
)
