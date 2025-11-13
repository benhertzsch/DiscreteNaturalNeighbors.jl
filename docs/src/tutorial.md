```@meta
CurrentModule = DiscreteNaturalNeighbors
```

# Tutorial

This tutorial demonstrates the use of the DiscreteNaturalNeighbors.jl package for the interpolation of a function from an unstructured grid in $n \gtreq 2$ dimensions.

## 2D interpolation

We first look at a 2D interpolation problem and define some function that we want to interpolate


```@example tutorial1
using DiscreteNaturalNeighbors
println("this is running")
```



```@example tutorial1
using Plots, DiscreteNaturalNeighbors

function f2D(x,y)
    return 0.1*(y-1)^2 / ( 2 + abs(x)^2)^2 + exp(-0.2*(y-x^2)^2) + exp(-0.4*(y-x)^2) * sin(0.25*(0.2*x^2+y)^2) + 1.2*exp(-((x+2)^2)/0.2 - (y-3.5)^4)
end

xRange, yRange = -5:0.1:5, -5:0.1:5

fun_vals = [f2D(x, y) for x in xRange, y in yRange]
heatmap(xRange, yRange, transpose(fun_vals), clim=(minimum(fun_vals), maximum(fun_vals)))
```

Next we generate a random sample of points on the 2D space with corresponding function values.

```@example tutorial1
NN     = 1000
points = 10 .* rand(NN, 2) .- 5;
values = [f2D(row...) for row in eachrow(points)];

println("size of points: ", size(points))
println("size of values: ", size(values))
```

Using these random samples, we then interpolate the function on a regular grid (i.e. mesh) over `(xRange, yRange)` by calling `interpolate_dnn()`:

```@example tutorial1
interp = interpolate_dnn(points, values, (xRange, yRange))

p1 = heatmap(xRange, yRange, transpose(interp), clim=(minimum(fun_vals), maximum(fun_vals)))
p2 = heatmap(xRange, yRange, transpose(interp - fun_vals))
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

The interpolation becomes more accurate for larger point sets. As the 2D interpolation is highly efficient even for large data sizes, progress is not displayed by default. A progress bar is enabled by setting the keyword argument `verbose=true`.

```@example tutorial1
NN     = 10000
points = 10 .* rand(NN, 2) .- 5;
values = [f2D(row...) for row in eachrow(points)];

interp = interpolate_dnn(points, values, (xRange, yRange); verbose=true)

p1 = heatmap(xRange, yRange, transpose(interp), clim=(minimum(fun_vals), maximum(fun_vals)))
p2 = heatmap(xRange, yRange, transpose(interp - fun_vals))
plot(p1, p2, layout=(1, 2), size=(800, 300))
```

**A note on efficiency:** Due to the inner loop structure of the DNN algorithm (see theory page), the interpolation is generally efficient when the target grid resolution is of the same order as the typical spacing between the sampled points; that is, if the number of points is comparable to the number of grid points.

## 3D interpolation

3D interpolations over a mesh `(xRange, yRange, zRange)` are obtained analogously to the 2D case. Again, we define some function to be interpolated, create random data points and use these to obtain the interpolation.

```@example tutorial1
function f3D(x, y, z)
    return exp(-z^2 * (1 + (0.1*x^2+0.4y^2))) * (0.1*(y-1)^2 / ( 2 + abs(x)^2)^2 + exp(-0.2*(y-x^2)^2) + exp(-0.4*(y-x)^2) * sin(0.25*(0.2*x^2+y)^2) + 1.2*exp(-((x+2)^2)/0.2 - (y-3.5)^4))
end

xRange, yRange, zRange = -5:0.1:5, -5:0.1:5, -5:0.1:5
#fun_vals = [f3D(x, y, z) for x in xRange, y in yRange, z in zRange]

NN     = 10000000
points = 10 .* rand(NN, 3) .- 5;
values = [f3D(row...) for row in eachrow(points)];


interp = interpolate_dnn(points, values, (xRange, yRange, zRange));

# plot a slice of the interpolation for comparison with the original function
z_idx = 44
fun_vals = [f3D(x, y, zRange[z_idx]) for x in xRange, y in yRange]
p1 = heatmap(xRange, yRange, fun_vals, clim=(minimum(fun_vals), maximum(fun_vals)))
p2 = heatmap(xRange, yRange, interp[:, :, z_idx], clim=(minimum(fun_vals), maximum(fun_vals)))
p3 = heatmap(xRange, yRange, fun_vals - interp[:, :, z_idx])
plot(p1, p2, p3, layout=(1, 3), size=(1200, 300))
```

Progress is shown for the 3D interpolation by default, and may be disabled by setting `verbose=false`.