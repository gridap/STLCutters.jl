module STLCuttersAnalysis

using DrWatson

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers
using GridapEmbedded
using STLCutters

import Downloads

using IterativeSolvers: cg
using Preconditioners: AMGPreconditioner
using Preconditioners: SmoothedAggregation
using STLCutters: FACE_IN, FACE_OUT, FACE_CUT

import STLCutters: check_requisites

export analysisdir
export testdir
export tmpdir

export download_thing
export download_run_and_save
export run_and_save
export rotations_and_displacements

include("Paths.jl")
include("RunSave.jl")

end

