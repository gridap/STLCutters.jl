module MPITestsBody

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("all_tests.jl")

if ! MPI.Initialized()
  MPI.Init()
end

if MPI.Comm_size(MPI.COMM_WORLD) == 8
  with_mpi() do distribute
    all_tests(distribute,(2,2,2))
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 1
  with_mpi() do distribute
    all_tests(distribute,(1,1,1))
  end
else
  MPI.Abort(MPI.COMM_WORLD,0)
end

end #module
