using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

include( scriptsdir("generic","create_jobs_matrix.jl") )
include( scriptsdir("generic","create_jobs_dataset.jl") )

for hpc_id in (:gadi,:titani,:acuario)

  create_jobs_matrix(hpc_id)

  create_jobs_dataset(hpc_id)

end

