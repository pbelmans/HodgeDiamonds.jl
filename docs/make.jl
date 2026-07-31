using AbstractAlgebra
using Documenter
using HodgeDiamonds

DocMeta.setdocmeta!(
    HodgeDiamonds, :DocTestSetup, :(using AbstractAlgebra, HodgeDiamonds); recursive=true
)

makedocs(;
    modules=[HodgeDiamonds],
    authors="Pieter Belmans",
    sitename="HodgeDiamonds.jl",
    format=Documenter.HTML(;
        canonical="https://pbelmans.github.io/HodgeDiamonds.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Hodge diamonds" => "diamonds.md",
        "Hochschild homology" => "hochschild.md",
        "Constructions" => "constructions.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/pbelmans/HodgeDiamonds.jl", devbranch="main")
