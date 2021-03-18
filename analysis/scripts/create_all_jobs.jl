using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

include( scriptsdir("generic","create_jobs_matrix.jl") )
include( scriptsdir("generic","create_jobs_dataset.jl") )

subsets = Dict(
  :acuario => 1:1000,
  :titani => 1001:3000,
  :gadi => 3001:5052 )

function create_all_jobs(hpc::Symbol,subset=subsets[hpc];onlymatrix=false,kwargs...)
  create_jobs_matrix(hpc;kwargs...)
  onlymatrix || create_jobs_dataset(hpc,subset;kwargs...)
  nothing
end
# create_jobs_matrix(:gadi,memory=32,displace=true,nmaxs=112,poisson=true,solution_order=1)
# create_jobs_matrix(:gadi,memory=32,displace=false,poisson=true,solution_order=2)
