using DrWatson

@quickactivate "STLCutters"

using STLCutters
using STLCutters.Tests

include(testdir("data","thingi10k_quality_filter.jl")) # file_ids = Int[]

download_and_run(file_ids[2:3],Δx=1e-5,vtk=true)

