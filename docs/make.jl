using Documenter, STLCutters

makedocs(;
    modules=[STLCutters],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pmartorell/STLCutters.jl/blob/{commit}{path}#L{line}",
    sitename="STLCutters.jl",
    authors="Pere Antoni Martorell, Large Scale Scientific Computing",
    assets=String[],
)

# deploydocs(;
#     repo="github.com/pmartorell/STLCutters.jl",
# )
