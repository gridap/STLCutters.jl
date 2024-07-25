
include("../CutterTests.jl")
include("../Poisson.jl")

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
