"""
    interpolate_dnn(points, values, ranges; verbose=nothing)

Evaluate the discrete natural neighbor interpolation of `points` (size `(N, n)`)
with `values` (size `(N,)`) on a regular grid (mesh) defined by
`ranges` (list of `n` ranges). If keyword argument `verbose` is unspecified,
default to no output for 2D interpolation, and progress bar for higher dimensions;
change by setting to `true` or `false`.
"""
function interpolate_dnn(points, values, ranges; verbose=nothing)
    @assert (typeof(points) == Matrix{Float64} && length(size(points)) == 2 && size(points)[2] > 1) (
        "Enter points as Matrix{Float64} of shape (N, ndim) with ndim > 1."
    )

    @assert (typeof(values) == Vector{Float64} && size(points)[1] == length(values)) (
        "Enter values as Vector{Float64} of length N matching first dimension of points (N, ndim)."
    )

    @assert length(ranges) == size(points)[2] (
        "Enter tuple of ranges (x1_Range, x2_Range, ..., xn_Range) matching dimensionality of points."
    )

    if size(points)[2] == 2   # 2D case
        if isnothing(verbose); verbose=false; end
        !verbose ? interpolate_dnn_2D(points, values, ranges) : interpolate_dnn_2D_(points, values, ranges);
    elseif size(points)[2] == 3   # 3D case
        if isnothing(verbose); verbose=true; end
        !verbose ? interpolate_dnn_3D(points, values, ranges) : interpolate_dnn_3D_(points, values, ranges);
    else   # nD case
        if isnothing(verbose); verbose=true; end
        !verbose ? interpolate_dnn_ND(points, values, ranges) : interpolate_dnn_ND_(points, values, ranges);
    end
    
end



###############################################################################
## 2D routine
###############################################################################

function interpolate_dnn_2D(points, values, (xRange, yRange))

    kdtree = KDTree(permutedims(points)) 

    Nx, Ny = length(xRange), length(yRange)

    grid       = [[x, y] for x in xRange, y in yRange]
    grid_val   = zeros(Nx, Ny)
    grid_count = zeros(Nx, Ny);

    step_x, step_y = xRange[2]-xRange[1], yRange[2]-yRange[1]

    for i in 1:Nx, j in 1:Ny
        p = grid[i, j]
        nn = knn(kdtree, p, 1, true)
        node_idx, node_dist = nn[1][1], nn[2][1]
        node = points[node_idx, :]

        # walk over minimal number of grid cells that may fall within dist < dist_node circle
        num_x, num_y = Int(ceil(node_dist / step_x)) + 1, Int(ceil(node_dist / step_y)) + 1
        sub_xRange = clamp(i - num_x, 1, Nx):clamp(i + num_x, 1, Nx)
        sub_yRange = clamp(j - num_y, 1, Ny):clamp(j + num_y, 1, Ny)

        for this_i in sub_xRange, this_j in sub_yRange
            this_candidate = grid[this_i, this_j]
            if (this_i == i) & (this_j == j) 
                grid_val[this_i, this_j]   += values[node_idx]
                grid_count[this_i, this_j] += 1
            else
                this_dist = norm(p, this_candidate)
                if this_dist <= node_dist
                    grid_val[this_i, this_j]   += values[node_idx]
                    grid_count[this_i, this_j] += 1
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end


function interpolate_dnn_2D_(points, values, (xRange, yRange))

    kdtree = KDTree(permutedims(points)) 

    Nx, Ny = length(xRange), length(yRange)

    grid       = [[x, y] for x in xRange, y in yRange]
    grid_val   = zeros(Nx, Ny)
    grid_count = zeros(Nx, Ny);

    step_x, step_y = xRange[2]-xRange[1], yRange[2]-yRange[1]

    @showprogress for i in 1:Nx
        for j in 1:Ny
            p = grid[i, j]
            nn = knn(kdtree, p, 1, true)
            node_idx, node_dist = nn[1][1], nn[2][1]
            node = points[node_idx, :]

            # walk over minimal number of grid cells that may fall within dist < dist_node circle
            num_x, num_y = Int(ceil(node_dist / step_x)) + 1, Int(ceil(node_dist / step_y)) + 1
            sub_xRange = clamp(i - num_x, 1, Nx):clamp(i + num_x, 1, Nx)
            sub_yRange = clamp(j - num_y, 1, Ny):clamp(j + num_y, 1, Ny)

            for this_i in sub_xRange, this_j in sub_yRange
                this_candidate = grid[this_i, this_j]
                if (this_i == i) & (this_j == j) 
                    grid_val[this_i, this_j]   += values[node_idx]
                    grid_count[this_i, this_j] += 1
                else
                    this_dist = norm(p, this_candidate)
                    if this_dist <= node_dist
                        grid_val[this_i, this_j]   += values[node_idx]
                        grid_count[this_i, this_j] += 1
                    end
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end



###############################################################################
## 3D routine
###############################################################################


function interpolate_dnn_3D(points, values, (xRange, yRange, zRange))

    kdtree = KDTree(permutedims(points)) 

    Nx, Ny, Nz = length(xRange), length(yRange), length(zRange)

    grid       = [[x, y, z] for x in xRange, y in yRange, z in zRange]
    grid_val   = zeros(Nx, Ny, Nz)
    grid_count = zeros(Nx, Ny, Nz);

    step_x, step_y, step_z = xRange[2]-xRange[1], yRange[2]-yRange[1], zRange[2]-zRange[1]

    for i in 1:Nx,j in 1:Ny, k in 1:Nz
        p = grid[i, j, k]
        nn = knn(kdtree, p, 1, true)
        node_idx, node_dist = nn[1][1], nn[2][1]
        node = points[node_idx, :]

        # walk over minimal number of grid cells that may fall within dist < dist_node sphere
        num_x, num_y, num_z = Int(ceil(node_dist / step_x)) + 1, Int(ceil(node_dist / step_y)) + 1, Int(ceil(node_dist / step_z)) + 1
        sub_xRange = clamp(i - num_x, 1, Nx):clamp(i + num_x, 1, Nx)
        sub_yRange = clamp(j - num_y, 1, Ny):clamp(j + num_y, 1, Ny)
        sub_zRange = clamp(k - num_z, 1, Nz):clamp(k + num_z, 1, Nz)

        for this_i in sub_xRange, this_j in sub_yRange, this_k in sub_zRange
            this_candidate = grid[this_i, this_j, this_k]
            
            if (this_i == i) & (this_j == j) & (this_k == k) 
                grid_val[this_i, this_j, this_k]   += values[node_idx]
                grid_count[this_i, this_j, this_k] += 1
            else
                this_dist = norm(p, this_candidate)
                if this_dist <= node_dist
                    grid_val[this_i, this_j, this_k]   += values[node_idx]
                    grid_count[this_i, this_j, this_k] += 1
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end


function interpolate_dnn_3D_(points, values, (xRange, yRange, zRange))

    kdtree = KDTree(permutedims(points)) 

    Nx, Ny, Nz = length(xRange), length(yRange), length(zRange)

    grid       = [[x, y, z] for x in xRange, y in yRange, z in zRange]
    grid_val   = zeros(Nx, Ny, Nz)
    grid_count = zeros(Nx, Ny, Nz);

    step_x, step_y, step_z = xRange[2]-xRange[1], yRange[2]-yRange[1], zRange[2]-zRange[1]

    @showprogress for i in 1:Nx
        for j in 1:Ny, k in 1:Nz
            p = grid[i, j, k]
            nn = knn(kdtree, p, 1, true)
            node_idx, node_dist = nn[1][1], nn[2][1]
            node = points[node_idx, :]

            # walk over minimal number of grid cells that may fall within dist < dist_node sphere
            num_x, num_y, num_z = Int(ceil(node_dist / step_x)) + 1, Int(ceil(node_dist / step_y)) + 1, Int(ceil(node_dist / step_z)) + 1
            sub_xRange = clamp(i - num_x, 1, Nx):clamp(i + num_x, 1, Nx)
            sub_yRange = clamp(j - num_y, 1, Ny):clamp(j + num_y, 1, Ny)
            sub_zRange = clamp(k - num_z, 1, Nz):clamp(k + num_z, 1, Nz)

            for this_i in sub_xRange, this_j in sub_yRange, this_k in sub_zRange
                this_candidate = grid[this_i, this_j, this_k]
                
                if (this_i == i) & (this_j == j) & (this_k == k) 
                    grid_val[this_i, this_j, this_k]   += values[node_idx]
                    grid_count[this_i, this_j, this_k] += 1
                else
                    this_dist = norm(p, this_candidate)
                    if this_dist <= node_dist
                        grid_val[this_i, this_j, this_k]   += values[node_idx]
                        grid_count[this_i, this_j, this_k] += 1
                    end
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end


###############################################################################
## ND routines
###############################################################################

function interpolate_dnn_ND(points, values, ranges)
    kdtree = KDTree(permutedims(points))

    ndim = size(points)[2]
    N = [length(range) for range in ranges]

    grid       = [collect(p) for p in Iterators.product(ranges...)]
    grid_val   = zeros(N...)
    grid_count = zeros(N...);

    step = [range[2]-range[1] for range in ranges]

    for grid_idx in CartesianIndices(tuple(N...))
        p  = grid[grid_idx]
        nn = knn(kdtree, p, 1, true)
        node_idx, node_dist = nn[1][1], nn[2][1]
        node = points[node_idx, :]

        num_steps = [Int(ceil(node_dist / step[d])) + 1 for d in 1:ndim]
        sub_ranges = [clamp(grid_idx[d] - num_steps[d], 1, N[d]):clamp(grid_idx[d] + num_steps[d], 1, N[d]) for d in 1:ndim]
        for this_idx in Iterators.product(sub_ranges...)
            this_idx = CartesianIndex(this_idx)
            p_prime = grid[this_idx]

            if this_idx == grid_idx
                grid_val[this_idx]   += values[node_idx]
                grid_count[this_idx] += 1
            else
                this_dist = norm(p, p_prime)
                if this_dist <= node_dist
                    grid_val[this_idx]   += values[node_idx]
                    grid_count[this_idx] += 1
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end


function interpolate_dnn_ND_(points, values, ranges)
    kdtree = KDTree(permutedims(points))

    ndim = size(points)[2]
    N = [length(range) for range in ranges]

    grid       = [collect(p) for p in Iterators.product(ranges...)]
    grid_val   = zeros(N...)
    grid_count = zeros(N...);

    step = [range[2]-range[1] for range in ranges]

    @showprogress for grid_idx in CartesianIndices(tuple(N...))
        p  = grid[grid_idx]
        nn = knn(kdtree, p, 1, true)
        node_idx, node_dist = nn[1][1], nn[2][1]
        node = points[node_idx, :]

        num_steps = [Int(ceil(node_dist / step[d])) + 1 for d in 1:ndim]
        sub_ranges = [clamp(grid_idx[d] - num_steps[d], 1, N[d]):clamp(grid_idx[d] + num_steps[d], 1, N[d]) for d in 1:ndim]
        for this_idx in Iterators.product(sub_ranges...)
            this_idx = CartesianIndex(this_idx)
            p_prime = grid[this_idx]

            if this_idx == grid_idx
                grid_val[this_idx]   += values[node_idx]
                grid_count[this_idx] += 1
            else
                this_dist = norm(p, p_prime)
                if this_dist <= node_dist
                    grid_val[this_idx]   += values[node_idx]
                    grid_count[this_idx] += 1
                end
            end
        end
    end

    interp = grid_val ./ grid_count
    return interp
end