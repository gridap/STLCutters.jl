using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

include( scriptsdir("generic","create_jobs_matrix.jl") )
include( scriptsdir("generic","create_jobs_dataset.jl") )

hpcs = (:gadi,:titani,:acuario)
subsets = [ 3001:5052, 1001:3000, 1:1000 ]

for (i,hpc_id) in enumerate(hpcs)

  create_jobs_matrix(hpc_id)

  create_jobs_dataset(hpc_id,subsets[i])

end

