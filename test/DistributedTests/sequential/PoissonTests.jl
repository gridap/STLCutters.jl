module PoissonTestsSeq
using PartitionedArrays
include("../Poisson.jl")
main = DistributedPoissonTests.main

with_debug() do distribute

  # Degenerated Case
  main(distribute,np=(1,1,1),nc=(4,4,4))

  # Distributed without propagation
  main(distribute,np=(2,2,2),nc=(4,4,4))
  main(distribute,np=(2,2,2),nc=(8,8,8))

end
end
