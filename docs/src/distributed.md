# Distributed

When dealing with large-scale problems, this package can be accelerate through two types of parallelization. The first one is multi-threading, which uses [Julia Threads](https://docs.julialang.org/en/v1/base/multi-threading/) for shared memory parallelization (e.g., `julia -t 4`). This method, adds some speed-up. However, it is only efficient for a reduced number of threads.

The second one is a distributed memory computing. For such parallelization, we use MPI (`mpiexec -np 4 julia input.jl`) though [`PartitionedArrays`](https://www.francescverdugo.com/PartitionedArrays.jl/stable) and [`GridapEmbedded`](https://gridap.github.io/GridapDistributed.jl/dev/). With MPI expect to efficiently compute large-scale problems, up to thousands of cores.


# Distributed usage

The STLCutters package has both multi-threading and distributed memory computing implementations. Multi-threading is straightforward for the user. However, setting distributed memory computations is a bit more involved at a driver level.

Here, we provide an example of a distributed computation with STLCutters.


!!! note
    Add example of distributed 
    We can also add a tutorial for the serial
    Add some explanation of the algorithm, cite PhD thesis
    Comment GridapEMbedded since 0.9.2



!!! warning
    Even though the distributed algorithms are proven to be efficient for large-scale weak scaling tests [Martorell, 2024]. The performance of this implementation is not tested.






