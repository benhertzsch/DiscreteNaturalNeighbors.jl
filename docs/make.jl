using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DiscreteNaturalNeighbors
using Documenter

DocMeta.setdocmeta!(DiscreteNaturalNeighbors, :DocTestSetup, :(using DiscreteNaturalNeighbors); recursive=true)

makedocs(;
    modules=[DiscreteNaturalNeighbors],
    authors="Benjamin Hertzsch <benjamin.hertzsch@ed.ac.uk> and contributors",
    sitename="DiscreteNaturalNeighbors.jl",
    format=Documenter.HTML(;
        canonical="https://benhertzsch.github.io/DiscreteNaturalNeighbors.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/benhertzsch/DiscreteNaturalNeighbors.jl",
    devbranch="main",
)
