using STLCutters
using STLCutters.Tests
using STLCutters.Tests: scriptsdir

include( scriptsdir("generic","create_jobs_matrix.jl") )
include( scriptsdir("generic","create_jobs_dataset.jl") )

subsets = Dict(
  :acuario => 1:1,
  :titani => 1:2924,
  :gadi => 2925:4963 )

function create_all_jobs(hpc::Symbol,subset=subsets[hpc];onlymatrix=false,kwargs...)
  create_jobs_matrix(hpc;kwargs...)
  onlymatrix || create_jobs_dataset(hpc,subset;kwargs...)
  nothing
end
# create_jobs_matrix(:gadi,memory=32,displace=true,nmaxs=112,poisson=true,solution_order=1)
# create_jobs_matrix(:gadi,memory=32,displace=false,poisson=true,solution_order=2)

# create_jobs_matrix(:gadi,walltime="1:00:00",memory=64,ncpus=16,displace=false,poisson=true,solution_order=2,solver=:amg)
