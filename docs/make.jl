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
    canonical="https://pbelmans.ncag.info/HodgeDiamonds.jl",
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
  # the test suite doctests the package root, which covers these pages too
  doctest=false,
)

deploydocs(; repo="github.com/pbelmans/HodgeDiamonds.jl", devbranch="main")
