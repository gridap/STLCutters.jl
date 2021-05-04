
include("CreateJobsMatrix.jl")

CreateJobsMatrix.main(:gadi,memory=32,ncpus=8,nmaxs=112,poisson=false)

CreateJobsMatrix.main(:gadi,memory=32,ncpus=8,walltime="5:00:00",
  displace=false,poisson=false)

CreateJobsMatrix.main(:gadi,memory=32,ncpus=8,
  nmaxs=112,poisson=true,solution_order=1,solver=:direct)

CreateJobsMatrix.main(:gadi,memory=32,ncpus=8,walltime="5:00:00",
  displace=false,poisson=true,solution_order=2,solver=:amg)

CreateJobsMatrix.main(:gadi,memory=16,ncpus=1,
  displace=false,poisson=true,solution_order=2,solver=:amg,nruns=5)

include("CreateJobsDataset.jl")

CreateJobsDataset.main(:gadi,2925:4963,ncpus=8,memory=32,chunk_size=100)

CreateJobsDataset.main(:titani,1:2924,ncpus=12,memory=120,chunk_size=500)
