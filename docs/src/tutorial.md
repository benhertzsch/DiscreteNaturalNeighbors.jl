```@meta
CurrentModule = DiscreteNaturalNeighbors
```

# Tutorial

This tutorial demonstrates the use of the DiscreteNaturalNeighbors.jl package for the interpolation of a function from an unstructured grid in $n \gtreq 2$ dimensions.

## 2D interpolation

We first look at a 2D interpolation problem and define some function that we want to interpolate

```@example tutorial1
using Plots, DiscreteNaturalNeighbors

function f2D(x,y)
    return 0.1*(y-1)^2 / ( 2 + abs(x)^2)^2 + exp(-0.2*(y-x^2)^2) + exp(-0.4*(y-x)^2) * sin(0.25*(0.2*x^2+y)^2) + 1.2*exp(-((x+2)^2)/0.2 - (y-3.5)^4)
end

xRange, yRange = -5:0.1:5, -5:0.1:5

fun_vals = [f2D(x, y) for x in xRange, y in yRange]
heatmap(xRange, yRange, transpose(fun_vals), clim=(minimum(fun_vals), maximum(fun_vals)))
```
