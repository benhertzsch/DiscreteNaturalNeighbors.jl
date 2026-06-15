using DiscreteNaturalNeighbors
using Test
using JLD2

@testset "DiscreteNaturalNeighbors.jl" begin
    
    # --- test 2D interpolation routines --- #
    println("running 2D interpolation tests")

    dict = load("data/data_2D.jld2")
    points, values = dict["points"], dict["values"]

    xRange, yRange = -5:0.05:5, -5:0.05:5
    interp   = interpolate_dnn(points, values, (xRange, yRange))
    @test size(interp) == (201, 201)
    @test interp[107, 48] ≈ 0.5557555333270575
    @test interp[133, 21] ≈ 0.1227873463298795

    interp_a = DiscreteNaturalNeighbors.interpolate_dnn_2D(points, values, (xRange, yRange));
    interp_b = DiscreteNaturalNeighbors.interpolate_dnn_ND(points, values, (xRange, yRange));
    @test size(interp_a) == (201, 201) && size(interp_b) == (201, 201)
    @test all(interp_a .≈ interp)
    @test all(interp_b .≈ interp)

    interp_a2 = DiscreteNaturalNeighbors.interpolate_dnn_2D_(points, values, (xRange, yRange));
    interp_b2 = DiscreteNaturalNeighbors.interpolate_dnn_ND_(points, values, (xRange, yRange));
    @test size(interp_a2) == (201, 201) && size(interp_b2) == (201, 201)
    @test all(interp_a2 .≈ interp)
    @test all(interp_b2 .≈ interp)


    # --- test 3D interpolation routines --- #
    println("running 3D interpolation tests")
    dict = load("data/data_3D.jld2")
    points, values = dict["points"], dict["values"]

    xRange, yRange, zRange = -5:0.1:5, -5:0.1:5, -5:0.1:5
    interp = interpolate_dnn(points, values, (xRange, yRange, zRange));
    @test size(interp) == (101, 101, 101)
    @test interp[64, 37, 51] ≈ 0.248817538389923

    interp_a = DiscreteNaturalNeighbors.interpolate_dnn_3D(points, values, (xRange, yRange, zRange));
    interp_b = DiscreteNaturalNeighbors.interpolate_dnn_ND(points, values, (xRange, yRange, zRange));
    @test size(interp_a) == (101, 101, 101) && size(interp_b) == (101, 101, 101)
    @test all(interp_a .≈ interp)
    @test all(interp_b .≈ interp)

end
