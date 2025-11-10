using DiscreteNaturalNeighbors
using Test

@testset "DiscreteNaturalNeighbors.jl" begin
    
    using JLD2

    # --- test 2D interpolation routines --- #
    dict = load("data/data_2D.jld2")
    points, values = dict["points"], dict["values"]

    xRange, yRange = -5:0.05:5, -5:0.05:5
    interp_a = interpolate2D(points, values, (xRange, yRange));
    interp_b = interpolateND(points, values, (xRange, yRange));

    @test size(interp_a) == (201, 201)
    @test size(interp_b) == (201, 201)
    @test interp_a[107, 48] ≈ 0.5589506815196893
    @test all(interp_a .≈ interp_b)


end
