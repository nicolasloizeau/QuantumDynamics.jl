
using QuantumDynamics
using Documenter

DocMeta.setdocmeta!(QuantumDynamics, :DocTestSetup, :(using QuantumDynamics); recursive=true)

makedocs(
    modules=[QuantumDynamics],
    authors="Nicolas Loizeau",
    sitename="QuantumDynamics.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://nicolasloizeau.github.io/QuantumDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Heisenberg time evolution" => "heisenberg.md",
    ],
)

deploydocs(
    repo = "github.com/nicolasloizeau/QuantumDynamics.jl.git",
    devbranch = "main",
    branch = "gh-pages",
)
