module Tests

using DrWatson

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers
using STLCutters

import HTTP

using STLCutters: FACE_IN, FACE_OUT, FACE_CUT

export analysisdir
export testdir
export tmpdir

export download_thing
export download_run_and_save
export run_and_save

include("Paths.jl")
include("RunSave.jl")

end

