module MPITestsBody

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("../CutterTests.jl")
include("../Poisson.jl")

if ! MPI.Initialized()
  MPI.Init()
end

function all_tests(distribute,parts)
  ranks = distribute(LinearIndices((prod(parts),)))

  t = PArrays.PTimer(ranks,verbose=true)
  PArrays.tic!(t)

  DistributedCutterTests.main(distribute,np=parts,nc=(4,4,4))
  DistributedCutterTests.main(distribute,np=parts,nc=(8,8,8))
  DistributedCutterTests.main(distribute,np=parts,nc=(4,4,4),simplex=true)
  DistributedCutterTests.main(distribute,np=parts,nc=(8,8,8),simplex=true)
  PArrays.toc!(t,"MPIDistributedCutter")

  DistributedPoissonTests.main(distribute,np=parts,nc=(4,4,4))
  DistributedPoissonTests.main(distribute,np=parts,nc=(8,8,8))
  DistributedPoissonTests.main(distribute,np=parts,nc=(8,8,8),geoname="Bunny-LowPoly")
  PArrays.toc!(t,"MPIDistributedPoisson")

  display(t)
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
