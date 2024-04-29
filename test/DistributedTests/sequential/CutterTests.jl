module CutterTestsSeq
using PartitionedArrays
include("../CutterTests.jl")
main = DistributedCutterTests.main

with_debug() do distribute

  # Degenerated Case
  main(distribute,np=(1,1,1),nc=(4,4,4))
  main(distribute,np=(1,1,1),nc=(4,4,4),simplex=true)

  # Distributed without propagation
  main(distribute,np=(2,2,2),nc=(4,4,4))
  main(distribute,np=(2,2,2),nc=(8,8,8))
  main(distribute,np=(2,2,2),nc=(4,4,4),simplex=true)
  main(distribute,np=(2,2,2),nc=(8,8,8),simplex=true)

  # Distributed with propagation
  main(distribute,np=(3,3,3),nc=(9,9,9))
  main(distribute,np=(3,3,3),nc=(9,9,9),simplex=true)

  # Untouched subdomain interiors
  main(distribute,np=(3,3,3),nc=(9,9,9),δ=0.5)
  main(distribute,np=(3,3,3),nc=(9,9,9),δ=0.5,simplex=true)

  # Embedded boundary on subdomain iterface
  main(distribute,np=(6,6,6),nc=(12,12,12),δ=0.5,simplex=true)

end
end
