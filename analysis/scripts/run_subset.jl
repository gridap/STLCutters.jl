using DrWatson

@quickactivate "STLCutters"

using STLCutters
using STLCutters.Tests

include(testdir("data","thingi10k_quality_filter.jl")) # file_ids = Int[]

#things=[37841,37881,39295,44395]
things=[96457,550964,293137]
#things = file_ids[2:3]

download_run_and_save(things,vtk=true)

