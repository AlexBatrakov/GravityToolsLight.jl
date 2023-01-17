abstract type General2DGrid end

struct Refinement2DGrid <: General2DGrid
    value::Dict{Symbol,Matrix{Float64}}
    status::Matrix{Int64}
    ref_level::Matrix{Int64}
    params::Dict{Symbol,Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    N_x::Int64
    N_y::Int64
    function Refinement2DGrid(value::Dict, status::Matrix{Int64}, ref_level::Matrix{Int64}, params::Dict, x::Vector{Float64}, y::Vector{Float64})
        return new(value, status, ref_level, params, x, y, length(x), length(y))
    end
end

function Refinement2DGrid(value::Dict, x::Vector{Float64}, y::Vector{Float64})
    params = Dict{Symbol,Float64}()
    status = fill(-1,length(x),length(y))
    ref_level = fill(0,length(x),length(y))
    return Refinement2DGrid(value, status, ref_level, params, x, y)
end

function precalculate_2DGrid(grid::Refinement2DGrid, target_function, params_function!)
    for i in 1:grid.N_x, j in 1:grid.N_y
        target_keys, target_values = target_function(grid.x[i], grid.y[j])
        for (i_key, key) in enumerate(target_keys)
            if haskey(grid.value, key)
                grid.value[key][i, j] = target_values[i_key]
            else
                grid.value[key] = fill(-1, grid.N_x, grid.N_y)
            end
        end
    end
    grid.status .= 1
    params_function!(grid)
    return grid
end

function parallel_precalculate_2DGrid(grid::Refinement2DGrid, target_function, params_function!)


    np = nprocs()  # determine the number of processes available

    target_keys, target_values = target_function(grid.x[1], grid.y[1], only_keys = true)
    for (i_key, key) in enumerate(target_keys)
        grid.value[key] = fill(-1, grid.N_x, grid.N_y)
    end
    i = 1
    j = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    function nextidx()
        idx = (i,j)
        if i < grid.N_x
            i += 1
        else 
            i = 1
            j += 1
        end
        return idx
    end
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx[2] > grid.N_y
                            break
                        end
                        println("p = $p, idx = $idx")
                        target_keys, target_values = remotecall_fetch(target_function, p, grid.x[idx[1]], grid.y[idx[2]])
                        for (i_key, key) in enumerate(target_keys)
                            grid.value[key][idx[1], idx[2]] = target_values[i_key]
                        end
                    end
                end
            end
        end
    end

    grid.status .= 1
    params_function!(grid)
    return grid
end