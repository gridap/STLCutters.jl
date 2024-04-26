module CutterTestsSeq
using PartitionedArrays
include("../CutterTests.jl")
main = DistributedCutterTests.main

with_debug() do distribute

  # Degenerated Case
  main(distribute,np=(1,1,1),nc=(4,4,4))
  # main(distribute,np=(1,1,1),nc=(4,4,4),simplex=true,verbose=true,vtk=true)

  # Distributed without propagation
  main(distribute,np=(2,2,2),nc=(4,4,4))
  main(distribute,np=(2,2,2),nc=(8,8,8))

  # Distributed with propagation
  main(distribute,np=(3,3,3),nc=(9,9,9))

  # Untouched subdomain interiors
  main(distribute,np=(3,3,3),nc=(9,9,9),δ=0.5)

  # Embedded boundary on subdomain iterface
  main(distribute,np=(6,6,6),nc=(12,12,12),δ=0.5)

end
end
