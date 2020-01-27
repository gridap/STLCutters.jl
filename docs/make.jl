using Documenter, STLCutter

makedocs(;
    modules=[STLCutter],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pmartorell/STLCutter.jl/blob/{commit}{path}#L{line}",
    sitename="STLCutter.jl",
    authors="Pere Antoni Martorell, Large Scale Scientific Computing",
    assets=String[],
)

deploydocs(;
    repo="github.com/pmartorell/STLCutter.jl",
)
