# STLCutter.jl

Welcome to the documentation for STLCutters.


!!! note

    These documentation pages are under construction.


## What

This package provides the tools for solving partial differential equations (PDEs) on domains bounded by [STL](https://en.wikipedia.org/wiki/STL_(file_format)) surfaces. STLCutters is an extension of [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl). GridapEmbedded is a package for solving 
PDEs on embedded domains, e.g., using embedded, unfitted or immersed finite element methods. Both, GridapEmbedded and STLCutters build on top of [Gridap](https://github.com/gridap/GridapEmbedded.jl).

In the following work, you can find the research related with STLCutters:

> Santiago Badia, Pere A. Martorell, Francesc Verdugo. "Geometrical discretisations for unfitted finite elements on explicit boundary representations." Journal of Computational Physics 460 (2022): 111162. doi: [10.1016/j.jcp.2022.111162](https://doi.org/10.1016/j.jcp.2022.111162)

## Why

The simulation of industrial and scientific problems involves solving PDEs on complex geometries. The generation of body-fitted unstructured meshes requires extensive human intervention. Additionaly, mesh partitioners are inherently serial and represent a botleneck (or a limitation) in the parallelization. Embedded methods (e.g., [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl)) address this limitation by using simple (e.g., structured) meshes for the functional discretization.. However, these methods define the domain through implicit functions (i.e., level sets) which represents a significant limitation.

This package addresses explicit boundary representations with embedded methods. Specifically, we provide algorithmic tools for the discretizations of embedded methods on STL surfaces. The implementations are designed to be efficient in large-scale distributed memory environments.
