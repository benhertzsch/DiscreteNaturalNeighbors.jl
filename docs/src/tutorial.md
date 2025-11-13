```@meta
CurrentModule = DiscreteNaturalNeighbors
```

# Tutorial

This tutorial demonstrates the use of the DiscreteNaturalNeighbors.jl package for the interpolation of a function from an unstructured grid in $n \gtreq 2$ dimensions.

## 2D interpolation

We first look at a 2D interpolation problem and define some function that we want to interpolate

```@example tutorial1
using DiscreteNaturalNeighbors, Plots
println("this is running")
```