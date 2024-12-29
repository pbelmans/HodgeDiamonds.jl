using HodgeDiamonds
using Documenter

DocMeta.setdocmeta!(HodgeDiamonds, :DocTestSetup, :(using HodgeDiamonds); recursive=true)

makedocs(;
    modules=[HodgeDiamonds],
    authors="Pieter Belmans <pieterbelmans@gmail.com> and contributors",
    sitename="HodgeDiamonds.jl",
    format=Documenter.HTML(;
        canonical="https://pbelmans.github.io/HodgeDiamonds.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pbelmans/HodgeDiamonds.jl",
    devbranch="main",
)
