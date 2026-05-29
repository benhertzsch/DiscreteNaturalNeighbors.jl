```@meta
CurrentModule = DiscreteNaturalNeighbors
```

# Tutorial

This tutorial demonstrates the use of the DiscreteNaturalNeighbors.jl package for the interpolation of a function from an unstructured grid in ``n \geq 2`` dimensions.

## 2D interpolation

We first look at a 2D interpolation problem and define some function that we want to interpolate

```@example tutorial1
using CairoMakie, DiscreteNaturalNeighbors

f2D(x,y) = 0.1*(y-1)^2 / ( 2 + abs(x)^2)^2 + exp(-0.2*(y-x^2)^2) + exp(-0.4*(y-x)^2) * sin(0.25*(0.2*x^2+y)^2) + 1.2*exp(-((x+2)^2)/0.2 - (y-3.5)^4)

xRange, yRange = -5:0.1:5, -5:0.1:5

fun_vals = [f2D(x, y) for x in xRange, y in yRange]

fig = Figure()
ax1 = Axis(fig[1, 1], aspect = 1)
heatmap!(ax1, xRange, yRange, transpose(fun_vals), colorrange=(minimum(fun_vals), maximum(fun_vals)))
fig
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

fig = Figure()
ax1 = Axis(fig[1, 1], aspect = 1)
ax2 = Axis(fig[1, 2], aspect = 1)

heatmap!(ax1, xRange, yRange, transpose(interp), colorrange=(minimum(fun_vals), maximum(fun_vals)))
heatmap!(ax2, xRange, yRange, transpose(interp - fun_vals))

fig
```

The interpolation becomes more accurate for larger point sets. As the 2D interpolation is highly efficient even for large data sizes, progress is not displayed by default. A progress bar is enabled by setting the keyword argument `verbose=true`.

```@example tutorial1
NN     = 10000
points = 10 .* rand(NN, 2) .- 5;
values = [f2D(row...) for row in eachrow(points)];

interp = interpolate_dnn(points, values, (xRange, yRange); verbose=true)

fig = Figure()
ax1 = Axis(fig[1, 1], aspect = 1)
ax2 = Axis(fig[1, 2], aspect = 1)

heatmap!(ax1, xRange, yRange, transpose(interp), colorrange=(minimum(fun_vals), maximum(fun_vals)))
heatmap!(ax2, xRange, yRange, transpose(interp - fun_vals))

fig
```

**A note on efficiency:** Due to the inner loop structure of the DNN algorithm (see theory page), the interpolation is generally efficient when the target grid resolution is of the same order as the typical spacing between the sampled points; that is, if the number of points is comparable to the number of grid points.

## 3D interpolation

3D interpolations over a mesh `(xRange, yRange, zRange)` are obtained analogously to the 2D case. Again, we define a function to be interpolated, create random data points and use these to obtain the interpolation.

```@example tutorial1
f3D(x, y, z) = exp(-z^2 * (1 + (0.1*x^2+0.4y^2))) * (0.1*(y-1)^2 / ( 2 + abs(x)^2)^2 + exp(-0.2*(y-x^2)^2) + exp(-0.4*(y-x)^2) * sin(0.25*(0.2*x^2+y)^2) + 1.2*exp(-((x+2)^2)/0.2 - (y-3.5)^4))

xRange, yRange, zRange = -5:0.1:5, -5:0.1:5, -5:0.1:5
#fun_vals = [f3D(x, y, z) for x in xRange, y in yRange, z in zRange]

NN     = 10000000
points = 10 .* rand(NN, 3) .- 5;
values = [f3D(row...) for row in eachrow(points)];

interp = interpolate_dnn(points, values, (xRange, yRange, zRange));

# plot a slice of the interpolation for comparison with the original function
z_idx = 44
fun_vals = [f3D(x, y, zRange[z_idx]) for x in xRange, y in yRange]

fig = Figure()
ax1 = Axis(fig[1, 1], aspect = 1)
ax2 = Axis(fig[1, 2], aspect = 1)
ax3 = Axis(fig[1, 3], aspect = 1)

heatmap!(ax1, xRange, yRange, fun_vals, colorrange=(minimum(fun_vals), maximum(fun_vals)))
heatmap!(ax2, xRange, yRange, interp[:, :, z_idx], colorrange=(minimum(fun_vals), maximum(fun_vals)))
heatmap!(ax3, xRange, yRange, fun_vals - interp[:, :, z_idx])

fig
```

Progress is shown for the 3D interpolation by default, and may be disabled by setting `verbose=false`.