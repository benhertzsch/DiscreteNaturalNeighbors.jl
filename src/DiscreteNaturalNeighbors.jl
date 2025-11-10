module DiscreteNaturalNeighbors

using LinearAlgebra, NearestNeighbors, ProgressMeter

include("interpolation.jl")

export interpolate2D, interpolate3D, interpolateND

end
