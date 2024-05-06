using Revise
using Documenter
using Literate
using STLCutters

makedocs(;
    modules=[STLCutters],
    format=Documenter.HTML(),
    pages=[
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "Distributed" => "distributed.md",
    ],
    sitename="STLCutters.jl",
)

deploydocs(
  repo = "github.com/gridap/STLCutters.jl.git",
)
