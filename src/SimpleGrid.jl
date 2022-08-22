struct SimpleGrid
    value::Dict{Symbol,Matrix{Float64}}
    status::Matrix{Int64}
    ref_level::Matrix{Int64}
    params::Dict{Symbol,Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    N_x::Int64
    N_y::Int64
    function SimpleGrid(value::Dict, status::Matrix{Int64}, ref_level::Matrix{Int64}, params::Dict, x::Vector{Float64}, y::Vector{Float64})
        return new(value, status, ref_level, params, x, y, length(x), length(y))
    end
end

function SimpleGrid(value::Dict, x::Vector{Float64}, y::Vector{Float64})
    params = Dict{Symbol,Float64}()
    status = fill(-1,length(x),length(y))
    ref_level = fill(0,length(x),length(y))
    return SimpleGrid(value, status, ref_level, params, x, y)
end

function precalculate_Grid(grid::SimpleGrid, get_value_func, calc_params!)
    sample_keys, sample_values = get_value_func(grid.x[1], grid.y[1])

    for key in sample_keys
        grid.value[key] = fill(-1, grid.N_x, grid.N_y)
    end
    for i in 1:grid.N_x, j in 1:grid.N_y
        calc_keys, calc_values = get_value_func(grid.x[i], grid.y[j])
        for (i_key, key) in enumerate(calc_keys)
            grid.value[key][i, j] = calc_values[i_key]
        end
    end
    grid.status .= 1
    calc_params!(grid)
    return grid
end

function refine_1Darray(x::Vector{Float64})
    x_refined = Vector{Float64}(undef, length(x)*2-1)
    for i in 1:length(x)-1
        x_refined[2*i-1] = x[i]
        x_refined[2*i] = 0.5*(x[i+1] + x[i])
    end
    x_refined[end] = x[end]
    return x_refined
end

function refine_2Darray(arr::Matrix{T}) where {T}
    arr_refined = fill(-one(T), 2 .* size(arr) .- 1)
    for i in 1:size(arr)[1], j in 1:size(arr)[2]
        arr_refined[2*i-1,2*j-1] = arr[i,j]
    end
    return arr_refined::Matrix{T}
end

function refine_Dict_of_2DArrays(dict::Dict{Symbol,Matrix{T}}) where {T}
    dict_refined = Dict{Symbol,Matrix{T}}()
    for (key, value) in dict
        dict_refined[key] = refine_2Darray(dict[key])
    end
    return dict_refined::Dict{Symbol,Matrix{T}}
end

function refine_Grid(grid::SimpleGrid, get_value_func, cell_selector, calc_params!)

    grid_refined = SimpleGrid(refine_Dict_of_2DArrays(grid.value), refine_2Darray(grid.status), refine_2Darray(grid.ref_level), grid.params, refine_1Darray(grid.x), refine_1Darray(grid.y))
    println("\nrefinement from ($(grid.N_x), $(grid.N_y)) to ($(grid_refined.N_x), $(grid_refined.N_y))")
    interp_counter = 0
    calc_counter = 0

    calc_status = fill(false, grid.N_x, grid.N_y)

#    for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1
#        if cell_selector(i_cell, j_cell, grid) == true
#            calculate_cell!(i_cell, j_cell, grid, grid_refined, get_value_func)
#            calc_status[i_cell,j_cell] = true
#            calc_counter += 1
#        else
#            interpolate_cell!(i_cell, j_cell, grid, grid_refined)
#            interp_counter += 1
#        end
#    end

    function recursive_refine_cell(i_cell, j_cell, calc_status, grid)
        if (i_cell ∉ 1:grid.N_x-1) || (j_cell ∉ 1:grid.N_y-1)
            return nothing
        end
        if (calc_status[i_cell,j_cell] == false) && cell_selector(i_cell, j_cell, grid) == true
            calculate_cell!(i_cell, j_cell, grid, grid_refined, get_value_func)
            calc_status[i_cell,j_cell] = true
            calc_counter += 1
            recursive_refine_cell(i_cell-1, j_cell, calc_status, grid)
            recursive_refine_cell(i_cell, j_cell-1, calc_status, grid)
            recursive_refine_cell(i_cell+1, j_cell, calc_status, grid)
            recursive_refine_cell(i_cell, j_cell+1, calc_status, grid)
        end
        return nothing
    end

#    changes_key = true
#    while changes_key
#        changes_key = false
        for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1
            recursive_refine_cell(i_cell, j_cell, calc_status, grid)
        end
#    end

    for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1
        if (calc_status[i_cell,j_cell] == false)
            interpolate_cell!(i_cell, j_cell, grid, grid_refined)
            interp_counter += 1
        end
    end


    println("calculations made : $calc_counter, interpolations made : $interp_counter")
    calc_params!(grid)
    return grid_refined
end

function interpolate_cell!(i::Int64, j::Int64, grid::SimpleGrid, grid_refined::SimpleGrid)
    i_ref = 2*i-1
    j_ref = 2*j-1
    
    if grid_refined.status[i_ref, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref, j_ref+1] = 0.5*(grid.value[key][i,j]+grid.value[key][i,j+1])
        end
        grid_refined.ref_level[i_ref, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i,j+1])
        grid_refined.status[i_ref, j_ref+1] = 0
    end
    if grid_refined.status[i_ref+1, j_ref] == -1
        for key in keys(grid_refined.value) 
            grid_refined.value[key][i_ref+1, j_ref] = 0.5*(grid.value[key][i,j]+grid.value[key][i+1,j])
        end
        grid_refined.ref_level[i_ref+1, j_ref] = min(grid.ref_level[i,j], grid.ref_level[i+1,j])
        grid_refined.status[i_ref+1, j_ref] = 0
    end
    if grid_refined.status[i_ref+1, j_ref+2] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+1, j_ref+2] = 0.5*(grid.value[key][i,j+1]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+2] = min(grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+2] = 0
    end
    if grid_refined.status[i_ref+2, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+2, j_ref+1] = 0.5*(grid.value[key][i+1,j]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+2, j_ref+1] = min(grid.ref_level[i+1,j], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+2, j_ref+1] = 0
    end
    if grid_refined.status[i_ref+1, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+1, j_ref+1] = 0.25*(grid.value[key][i,j]+grid.value[key][i+1,j]+grid.value[key][i,j+1]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i+1,j], grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+1] = 0
    end
    return grid_refined
end

function calculate_cell!(i_cell::Int64, j_cell::Int64, grid::SimpleGrid, grid_refined::SimpleGrid, get_value_func)
    i_cell_ref = 2*i_cell-1
    j_cell_ref = 2*j_cell-1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1

    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        if grid_refined.status[i_ref,j_ref] < 1
            is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1) 
            calc_keys, calc_values = get_value_func(grid_refined.x[i_ref], grid_refined.y[j_ref])
            for (i_key, key) in enumerate(calc_keys)
                grid_refined.value[key][i_ref, j_ref] = calc_values[i_key]
                if is_on_grid
                    grid.value[key][div(i_ref,2)+1,div(j_ref,2)+1] = calc_values[i_key]
                    grid.ref_level[div(i_ref,2)+1,div(j_ref,2)+1] = new_ref_level-1
                end
            end
            grid_refined.status[i_ref,j_ref] = 1
        end
        if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
            grid_refined.ref_level[i_ref, j_ref] = new_ref_level
        end
    end

    return grid_refined
end

grid_size_counter(;grid_init,ref_level) = (grid_init-1)*2^ref_level+1