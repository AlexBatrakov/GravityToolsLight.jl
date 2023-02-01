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

function simple_parallel_precalculate_2DGrid(grid::Refinement2DGrid, target_function, params_function!)

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

    n_steps = grid.N_x*grid.N_y
    p = Progress(n_steps)
    channel = RemoteChannel(()->Channel{Bool}(), 1)

    @sync begin
        @async while take!(channel)
            next!(p)
        end
        @sync for p=1:np
            if p != myid() || np == 1
                @async while true
                    idx = nextidx()
                    if idx[2] > grid.N_y
                        break
                    end
#                    println("myid = $(myid()), p = $p, idx = $idx")
                    target_keys, target_values = remotecall_fetch(target_function, p, grid.x[idx[1]], grid.y[idx[2]])
                    put!(channel, true)
                    for (i_key, key) in enumerate(target_keys)
                        grid.value[key][idx[1], idx[2]] = target_values[i_key]
                    end
                end
            end
        end
        put!(channel, false)
    end

    grid.status .= 1
    params_function!(grid)
    return grid
end

function pmap_parallel_precalculate_2DGrid(grid::Refinement2DGrid, target_function, params_function!)

    target_keys, target_values = target_function(grid.x[1], grid.y[1], only_keys = true)
    for (i_key, key) in enumerate(target_keys)
        grid.value[key] = fill(-1, grid.N_x, grid.N_y)
    end

    function separate_task(idx::Tuple{Int64,Int64})
        i, j = idx
        calc_keys, calc_values = target_function(grid.x[i], grid.y[j])
        for (i_key, key) in enumerate(calc_keys)
            grid.value[key][i, j] = calc_values[i_key]
        end
    end

    ans = pmap(args -> target_function(args...), [(x,y) for x in grid.x, y in grid.y])

    for i in 1:grid.N_x, j in 1:grid.N_y
        calc_keys, calc_values = ans[i,j]
        for (i_key, key) in enumerate(calc_keys)
            grid.value[key][i, j] = calc_values[i_key]
        end
    end

    grid.status .= 1
    params_function!(grid)
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

#--------------------------------------------------------------------------------------------------------------

function parallel_refine_2DGrid(grid::Refinement2DGrid, target_function, cell_selector, params_function!)

    np = nprocs()

    grid_refined = Refinement2DGrid(refine_Dict_of_2DArrays(grid.value), refine_2Darray(grid.status), refine_2Darray(grid.ref_level), grid.params, refine_1Darray(grid.x), refine_1Darray(grid.y))
    println("\nRefinement from ($(grid.N_x), $(grid.N_y)) to ($(grid_refined.N_x), $(grid_refined.N_y))")
    interp_counter = 0
    calc_counter = 0

    calc_status = fill(false, grid.N_x, grid.N_y)
    cell_selector_status = [true for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1]

    interations_counter = 0

    while true


    cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
    for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1
        if cell_selector(i_cell, j_cell, grid) && cell_selector_status[i_cell, j_cell]
            push!(cells_to_refine, (i_cell, j_cell))
        end
    end

    n_cells = length(cells_to_refine)

    if n_cells == 0
        break
    end
    interations_counter += 1
    counter_cells = 1

    nextidx() = (idx=counter_cells; counter_cells+=1; idx)

    p = Progress(n_cells)
    channel = RemoteChannel(()->Channel{Bool}(), 1)

    @sync begin
        @async while take!(channel)
            next!(p)
        end
        @sync begin
            for p=1:np
                if p != myid() || np == 1
                    @async begin
                        while true
                            idx = nextidx()
                            if idx > n_cells
                                break
                            end
                            i_cell, j_cell = cells_to_refine[idx]
                            calc_counter += calculate_cell!(p, i_cell, j_cell, grid, grid_refined, target_function)
                            cell_selector_status[i_cell, j_cell] = false
                            put!(channel, true)
                        end
                    end
                end
            end
        end
        put!(channel, false)
    end

    #interpolate the rest

    for i_cell in 1:grid.N_x-1, j_cell in 1:grid.N_y-1
        if (calc_status[i_cell,j_cell] == false)
            interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
        end
    end

    end

    println("interations = $interations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)
    return grid_refined
end

function interpolate_cell!(i::Int64, j::Int64, grid::Refinement2DGrid, grid_refined::Refinement2DGrid)
    i_ref = 2*i-1
    j_ref = 2*j-1
    inter_counter = 0
    if grid_refined.status[i_ref, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref, j_ref+1] = 0.5*(grid.value[key][i,j]+grid.value[key][i,j+1])
        end
        grid_refined.ref_level[i_ref, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i,j+1])
        grid_refined.status[i_ref, j_ref+1] = 0
        inter_counter += 1
    end
    if grid_refined.status[i_ref+1, j_ref] == -1
        for key in keys(grid_refined.value) 
            grid_refined.value[key][i_ref+1, j_ref] = 0.5*(grid.value[key][i,j]+grid.value[key][i+1,j])
        end
        grid_refined.ref_level[i_ref+1, j_ref] = min(grid.ref_level[i,j], grid.ref_level[i+1,j])
        grid_refined.status[i_ref+1, j_ref] = 0
        inter_counter += 1
    end
    if grid_refined.status[i_ref+1, j_ref+2] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+1, j_ref+2] = 0.5*(grid.value[key][i,j+1]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+2] = min(grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+2] = 0
        inter_counter += 1
    end
    if grid_refined.status[i_ref+2, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+2, j_ref+1] = 0.5*(grid.value[key][i+1,j]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+2, j_ref+1] = min(grid.ref_level[i+1,j], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+2, j_ref+1] = 0
        inter_counter
    end
    if grid_refined.status[i_ref+1, j_ref+1] == -1 
        for key in keys(grid_refined.value)
            grid_refined.value[key][i_ref+1, j_ref+1] = 0.25*(grid.value[key][i,j]+grid.value[key][i+1,j]+grid.value[key][i,j+1]+grid.value[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i+1,j], grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+1] = 0
        inter_counter
    end
    return inter_counter
end

function calculate_cell!(p, i_cell::Int64, j_cell::Int64, grid::Refinement2DGrid, grid_refined::Refinement2DGrid, target_function)
    i_cell_ref = 2*i_cell-1
    j_cell_ref = 2*j_cell-1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1
    calc_counter = 0
    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        if grid_refined.status[i_ref,j_ref] < 1
            is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1) 
            calc_keys, calc_values = remotecall_fetch(target_function, p, grid_refined.x[i_ref], grid_refined.y[j_ref])
            calc_counter +=1
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

    return calc_counter
end