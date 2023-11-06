#--------------------------------------------------------------------------------------------------------------
# Refinement units

abstract type AbstractRefinementUnit end

struct FullUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
end

FullUnit(name; min = -Inf, max = Inf) = FullUnit(name, min, max)

function Base.show(io::IO, ru::FullUnit)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Full refinement unit:")
    println(io, ' '^(indent+4), "Name of variable: ", ru.name)
    println(io, ' '^(indent+4), "Minimum value: ", ru.min)
    print(io,   ' '^(indent+4), "Maximum value: ", ru.max)
	return nothing
end

struct DiffUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diff::Float64
end

DiffUnit(name; min = -Inf, max = Inf, diff = diff) = DiffUnit(name, min, max, diff)

function Base.show(io::IO, ru::DiffUnit)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Difference refinement unit:")
    println(io, ' '^(indent+4), "Name of variable: ", ru.name)
    println(io, ' '^(indent+4), "Minimum value: ", ru.min)
    println(io, ' '^(indent+4), "Maximum value: ", ru.max)
    print(io,   ' '^(indent+4), "Maximal difference: ", ru.diff)
	return nothing
end

struct ContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    contours::Vector{Float64}
end

ContourUnit(name; min = -Inf, max = Inf, contours = contours) = ContourUnit(name, min, max, contours)

function Base.show(io::IO, ru::ContourUnit)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Contour refinement unit:")
    println(io, ' '^(indent+4), "Name of variable: ", ru.name)
    println(io, ' '^(indent+4), "Minimum value: ", ru.min)
    println(io, ' '^(indent+4), "Maximum value: ", ru.max)
    print(io,   ' '^(indent+4), "Contour levels: ", ru.contours)
	return nothing
end

struct DiffContourUnit <: AbstractRefinementUnit
    name::Symbol
    min::Float64
    max::Float64
    diff::Float64
    contours::Vector{Float64}
end

DiffContourUnit(name; min = -Inf, max = Inf, diff = diff, contours = contours) = DiffContourUnit(name, min, max, diff, contours)

function Base.show(io::IO, ru::DiffContourUnit)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Difference and contour refinement unit:")
    println(io, ' '^(indent+4), "Name of variable: ", ru.name)
    println(io, ' '^(indent+4), "Minimum value: ", ru.min)
    println(io, ' '^(indent+4), "Maximum value: ", ru.max)
    println(io, ' '^(indent+4), "Maximal difference: ", ru.diff)
    print(io,   ' '^(indent+4), "Contour levels: ", ru.contours)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Refinement settings

struct RefinementSettings{T}
    desired_refinement_level::Int64
    parallel::Bool
    units::T
end

function Base.show(io::IO, ref_sets::RefinementSettings{T}) where {T}
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "Grid Refinement settings:")
	println(io, ' '^(indent+4), "Desired refinement level: ", ref_sets.desired_refinement_level)
    println(io, ' '^(indent+4), "Parallel computation: ", ref_sets.parallel)
    for i in 1:length(ref_sets.units)-1
        println(IOContext(io, :indent => indent+4), ref_sets.units[i])
    end
    print(IOContext(io, :indent => indent+4), ref_sets.units[end])
	return nothing
end

function RefinementSettings(units...; desired_refinement_level::Int64, parallel::Bool)
    return RefinementSettings(desired_refinement_level, parallel, units)
end

#--------------------------------------------------------------------------------------------------------------
# Refinement 2D grid

abstract type General2DGrid end

struct AdaptiveRefinement2DGrid{T1 <: AbstractRangeRule, T2 <: AbstractRangeRule, T3} <: General2DGrid
    vars::Dict{Symbol,Matrix{Float64}}
    params::Dict{Symbol,Float64}
    min::Dict{Symbol,Float64}
    max::Dict{Symbol,Float64}
    x::RangeVariable{T1}
    y::RangeVariable{T2}
    ref_sets::RefinementSettings{T3}
    ref_level::Matrix{Int64}
    status::Matrix{Int64}
end

function AdaptiveRefinement2DGrid(x::RangeVariable{T1}, y::RangeVariable{T2}, ref_sets::T3) where {T1 <: AbstractRangeRule, T2 <: AbstractRangeRule, T3}
    vars = Dict{Symbol,Matrix{Float64}}()
    params = Dict{Symbol,Float64}()
    min = Dict{Symbol,Float64}()
    max = Dict{Symbol,Float64}()
    ref_level = [0 for i in 1:x.N, j in 1:y.N]
    status = [-1 for i in 1:x.N, j in 1:y.N]
    return AdaptiveRefinement2DGrid(vars, params, min, max, x, y, ref_sets, ref_level, status)
end

function Base.show(io::IO, grid::AdaptiveRefinement2DGrid)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "AdaptiveRefinement2DGrid:")
	println(io, ' '^(indent+4), "Variables: ", keys(grid.vars))
    println(io, ' '^(indent+4), "Parameters: ", grid.params)
    println(io, ' '^(indent+4), "Minimal values: ", grid.min)
    println(io, ' '^(indent+4), "Maximal values: ", grid.max)
    println(io, ' '^(indent+4), "X axis: ", grid.x)
    println(io, ' '^(indent+4), "Y axis: ", grid.y)
    print(IOContext(io, :indent => indent+4), grid.ref_sets)
	return nothing
end

#--------------------------------------------------------------------------------------------------------------
# Cell selectors for refinement units

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid)
    if !((0 < i_cell < grid.x.N) && (0 < j_cell < grid.y.N))
        return false
    end
    combined_case = false
    for ref_unit in grid.ref_sets.units
        combined_case = combined_case || cell_selector(i_cell, j_cell, grid, ref_unit)
    end
    return combined_case::Bool
end

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::FullUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    return min_case && max_case
end

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::DiffUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    diff_case = value_cell_max - value_cell_min > ref_unit.diff
    return diff_case && min_case && max_case
end
    
function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::ContourUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    value_min = grid.min[ref_unit.name]
    contour_case = any(value_cell_min .< value_min .+ ref_unit.contours .< value_cell_max)
    return contour_case && min_case && max_case
end

function cell_selector(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, ref_unit::DiffContourUnit)
    cell = @view grid.vars[ref_unit.name][i_cell:i_cell+1,j_cell:j_cell+1]
    value_cell_min = minimum(cell)
    value_cell_max = maximum(cell)
    min_case = ref_unit.min <= value_cell_min
    max_case = ref_unit.max >= value_cell_max
    value_min = grid.min[ref_unit.name]
    contour_case = any(value_cell_min .< value_min .+ ref_unit.contours .< value_cell_max)
    diff_case = (value_cell_max - value_cell_min > ref_unit.diff)
    return diff_case && contour_case && min_case && max_case
end

#--------------------------------------------------------------------------------------------------------------
# Supplementary refinement routines

function refine(var::RangeVariable)
    var_refined = typeof(var)(var.name, var.min, var.max, var.N*2-1, var.range_rule)
    return var_refined
end

#=
function refine_1Darray(x::Vector{Float64})
    x_refined = Vector{Float64}(undef, length(x)*2-1)
    for i in 1:length(x)-1
        x_refined[2*i-1] = x[i]
        x_refined[2*i] = 0.5*(x[i+1] + x[i])
    end
    x_refined[end] = x[end]
    return x_refined
end
=#

function refine(arr::Matrix{T}) where {T}
    arr_refined = fill(-one(T), 2 .* size(arr) .- 1)
    for i in 1:size(arr)[1], j in 1:size(arr)[2]
        arr_refined[2*i-1,2*j-1] = arr[i,j]
    end
    return arr_refined::Matrix{T}
end

function refine(dict::Dict{Symbol,Matrix{T}}) where {T}
    dict_refined = Dict{Symbol,Matrix{T}}()
    for (key, value) in dict
        dict_refined[key] = refine(dict[key])
    end
    return dict_refined::Dict{Symbol,Matrix{T}}
end

function refine(grid::AdaptiveRefinement2DGrid)
    vars_refined = refine(grid.vars)
    params_refined = copy(grid.params)
    min_refined = copy(grid.min)
    max_refined = copy(grid.max)
    x_refined = refine(grid.x)
    y_refined = refine(grid.y)
    ref_sets_refined = grid.ref_sets
    ref_level_refined = refine(grid.ref_level)
    status_refined = refine(grid.status)
    grid_refined = AdaptiveRefinement2DGrid(vars_refined, params_refined, min_refined, max_refined, x_refined, y_refined, ref_sets_refined, ref_level_refined, status_refined)
    return grid_refined
end

#--------------------------------------------------------------------------------------------------------------
# General routines

function calculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    precalculate_2DGrid!(grid, target_function, params_function!)
    for i in 1:grid.ref_sets.desired_refinement_level
        grid = refine_2DGrid(grid, target_function, params_function!)
    end
    return grid
end

function precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_precalculate_2DGrid!(grid, target_function, params_function!)
    else 
        return single_core_precalculate_2DGrid!(grid, target_function, params_function!)
    end
end

function refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    if grid.ref_sets.parallel == true
        return parallel_refine_2DGrid(grid, target_function, params_function!)
    else 
        return single_core_refine_2DGrid(grid, target_function, params_function!)
    end
end

#--------------------------------------------------------------------------------------------------------------
# Singe-core routines

function single_core_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    target_keys = Base.return_types(target_function, Tuple{Float64,Float64})[1].parameters[1]
    for key in target_keys
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end
    for i in 1:grid.x.N, j in 1:grid.y.N
        target_output = target_function(grid.x.values[i], grid.y.values[j])
        for (key, value) in pairs(target_output)
            grid.vars[key][i, j] = value
        end
    end
    for key in target_keys
        grid.min[key] = minimum(x->isnan(x) ? +Inf : x, grid.vars[key])
        grid.max[key] = maximum(x->isnan(x) ? -Inf : x, grid.vars[key])
    end
    grid.status .= 1
    params_function!(grid)
    return grid
end

function single_core_refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    grid_refined = refine(grid)

    println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
    interp_counter = 0
    calc_counter = 0
    iterations_counter = 0

    cell_selector_status = fill(false, grid.x.N-1, grid.y.N-1)

    while true
        cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
        for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
            if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                push!(cells_to_refine, (i_cell, j_cell))
                cell_selector_status[i_cell, j_cell] = true
            end
        end

        n_cells = length(cells_to_refine)

        if n_cells == 0
            break
        end
        iterations_counter += 1

        for cell in cells_to_refine
            i_cell, j_cell = cell
            calc_counter += calculate_cell!(i_cell, j_cell, grid, grid_refined, target_function)
            for cell in [(i_cell-1, j_cell), (i_cell, j_cell-1), (i_cell+1, j_cell), (i_cell, j_cell+1)]
                i_cell, j_cell = cell
                if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                    push!(cells_to_refine, cell)
                    cell_selector_status[i_cell, j_cell] = true
                end
            end
        end
    end

    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        if (cell_selector_status[i_cell,j_cell] == false)
            interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
        end
    end

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)
    return grid_refined
end

function calculate_cell!(i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, target_function)
    i_cell_ref = 2*i_cell-1
    j_cell_ref = 2*j_cell-1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1
    n_calc = 0

    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        if grid_refined.status[i_ref,j_ref] < 1
            is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1) 
            target_output = target_function(grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])
            n_calc += 1
            for (key, value) in pairs(target_output)
                grid_refined.vars[key][i_ref, j_ref] = value
                if !isnan(value)
                    grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
                    grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
                end
                if is_on_grid
                    grid.vars[key][div(i_ref,2)+1,div(j_ref,2)+1] = value
                    grid.ref_level[div(i_ref,2)+1,div(j_ref,2)+1] = new_ref_level-1
                    if !isnan(value)
                        grid.min[key] = value < grid.min[key] ? value : grid.min[key]
                        grid.max[key] = value > grid.max[key] ? value : grid.max[key]
                    end
                end
            end
            grid_refined.status[i_ref,j_ref] = 1
        end
        if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
            grid_refined.ref_level[i_ref, j_ref] = new_ref_level
        end
    end

    return n_calc
end

function interpolate_cell!(i::Int64, j::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid)
    i_ref = 2*i-1
    j_ref = 2*j-1
    n_inter = 0
    if grid_refined.status[i_ref, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref, j_ref+1] = 0.5*(grid.vars[key][i,j]+grid.vars[key][i,j+1])
        end
        grid_refined.ref_level[i_ref, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i,j+1])
        grid_refined.status[i_ref, j_ref+1] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+1, j_ref] == -1
        for key in keys(grid_refined.vars) 
            grid_refined.vars[key][i_ref+1, j_ref] = 0.5*(grid.vars[key][i,j]+grid.vars[key][i+1,j])
        end
        grid_refined.ref_level[i_ref+1, j_ref] = min(grid.ref_level[i,j], grid.ref_level[i+1,j])
        grid_refined.status[i_ref+1, j_ref] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+1, j_ref+2] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+1, j_ref+2] = 0.5*(grid.vars[key][i,j+1]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+2] = min(grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+2] = 0
        n_inter += 1
    end
    if grid_refined.status[i_ref+2, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+2, j_ref+1] = 0.5*(grid.vars[key][i+1,j]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+2, j_ref+1] = min(grid.ref_level[i+1,j], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+2, j_ref+1] = 0
        n_inter
    end
    if grid_refined.status[i_ref+1, j_ref+1] == -1 
        for key in keys(grid_refined.vars)
            grid_refined.vars[key][i_ref+1, j_ref+1] = 0.25*(grid.vars[key][i,j]+grid.vars[key][i+1,j]+grid.vars[key][i,j+1]+grid.vars[key][i+1,j+1])
        end
        grid_refined.ref_level[i_ref+1, j_ref+1] = min(grid.ref_level[i,j], grid.ref_level[i+1,j], grid.ref_level[i,j+1], grid.ref_level[i+1,j+1])
        grid_refined.status[i_ref+1, j_ref+1] = 0
        n_inter += 1
    end
    return n_inter
end


#--------------------------------------------------------------------------------------------------------------
# Parallel routines

function parallel_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    np = nprocs()  # determine the number of processes available
    target_keys = Base.return_types(target_function, Tuple{Float64,Float64})[1].parameters[1]
    for key in target_keys
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end

    i = 1
    j = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    function nextidx()
        idx = (i,j)
        if i < grid.x.N
            i += 1
        else 
            i = 1
            j += 1
        end
        return idx
    end

    n_steps = grid.x.N*grid.y.N
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
                    if idx[2] > grid.y.N
                        break
                    end
#                    println("myid = $(myid()), p = $p, idx = $idx")
                    target_output = remotecall_fetch(target_function, p, grid.x.values[idx[1]], grid.y.values[idx[2]])
                    put!(channel, true)
                    for (key, value) in pairs(target_output)
                        grid.vars[key][idx[1], idx[2]] = value
                    end
                end
            end
        end
        put!(channel, false)
    end

    for key in target_keys
        grid.min[key] = minimum(x->isnan(x) ? +Inf : x, grid.vars[key])
        grid.max[key] = maximum(x->isnan(x) ? -Inf : x, grid.vars[key])
    end
    grid.status .= 1
    params_function!(grid)
    return grid
end

#=
function pmap_parallel_precalculate_2DGrid!(grid::AdaptiveRefinement2DGrid, target_function, params_function!)

    target_keys, target_values = target_function(grid.x.values[1], grid.y.values[1], only_keys = true)
    for (i_key, key) in enumerate(target_keys)
        grid.vars[key] = fill(-1, grid.x.N, grid.y.N)
    end

    function separate_task(idx::Tuple{Int64,Int64})
        i, j = idx
        calc_keys, calc_values = target_function(grid.x.values[i], grid.y.values[j])
        for (i_key, key) in enumerate(calc_keys)
            grid.vars[key][i, j] = calc_values[i_key]
        end
    end

    ans = pmap(args -> target_function(args...), [(x,y) for x in grid.x, y in grid.y])

    for i in 1:grid.x.N, j in 1:grid.y.N
        calc_keys, calc_values = ans[i,j]
        for (i_key, key) in enumerate(calc_keys)
            grid.vars[key][i, j] = calc_values[i_key]
        end
    end

    grid.status .= 1
    params_function!(grid)
    return grid
end
=#


function parallel_refine_2DGrid(grid::AdaptiveRefinement2DGrid, target_function, params_function!)
    grid_refined = refine(grid)
    println("\nRefinement from ($(grid.x.N), $(grid.y.N)) to ($(grid_refined.x.N), $(grid_refined.y.N))")
    interp_counter = 0
    calc_counter = 0
    iterations_counter = 0
    np = nprocs()

    update_calc_counter(n_calc) = (calc_counter+=n_calc)

    cell_selector_status = fill(false, grid.x.N-1, grid.y.N-1)

    while true
        cells_to_refine = Vector{Tuple{Int64, Int64}}(undef, 0)
        for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
            if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                push!(cells_to_refine, (i_cell, j_cell))
                cell_selector_status[i_cell, j_cell] = true
            end
        end

        n_cells = length(cells_to_refine)

        if n_cells == 0
            break
        end
        iterations_counter += 1

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
                                n_calc = calculate_cell!(p, i_cell, j_cell, grid, grid_refined, target_function)
                                update_calc_counter(n_calc)
                                for cell in [(i_cell-1, j_cell), (i_cell, j_cell-1), (i_cell+1, j_cell), (i_cell, j_cell+1)]
                                    i_cell, j_cell = cell
                                    if cell_selector(i_cell, j_cell, grid) && !cell_selector_status[i_cell, j_cell]
                                        push!(cells_to_refine, cell)
                                        cell_selector_status[i_cell, j_cell] = true
                                    end
                                end
                                put!(channel, true)
                            end
                        end
                    end
                end
            end
            put!(channel, false)
        end
    end

    for i_cell in 1:grid.x.N-1, j_cell in 1:grid.y.N-1
        if (cell_selector_status[i_cell,j_cell] == false)
            interp_counter += interpolate_cell!(i_cell, j_cell, grid, grid_refined)
        end
    end

    println("iterations = $iterations_counter, calculations = $calc_counter, interpolations = $interp_counter")
    params_function!(grid)
    return grid_refined
end

function calculate_cell!(p, i_cell::Int64, j_cell::Int64, grid::AdaptiveRefinement2DGrid, grid_refined::AdaptiveRefinement2DGrid, target_function)
    i_cell_ref = 2*i_cell-1
    j_cell_ref = 2*j_cell-1
    new_ref_level = maximum(grid.ref_level[i_cell:i_cell+1,j_cell:j_cell+1]) + 1
    n_calc = 0

    for i_ref in i_cell_ref:i_cell_ref+2, j_ref in j_cell_ref:j_cell_ref+2
        if grid_refined.status[i_ref,j_ref] < 1
            is_on_grid = (mod(i_ref,2)*mod(j_ref,2) == 1) 
            target_output = remotecall_fetch(target_function, p, grid_refined.x.values[i_ref], grid_refined.y.values[j_ref])
            n_calc += 1
            for (key, value) in pairs(target_output)
                grid_refined.vars[key][i_ref, j_ref] = value
                if !isnan(value)
                    grid_refined.min[key] = value < grid_refined.min[key] ? value : grid_refined.min[key]
                    grid_refined.max[key] = value > grid_refined.max[key] ? value : grid_refined.max[key]
                end
                if is_on_grid
                    grid.vars[key][div(i_ref,2)+1,div(j_ref,2)+1] = value
                    grid.ref_level[div(i_ref,2)+1,div(j_ref,2)+1] = new_ref_level-1
                    if !isnan(value)
                        grid.min[key] = value < grid.min[key] ? value : grid.min[key]
                        grid.max[key] = value > grid.max[key] ? value : grid.max[key]
                    end
                end
            end
            grid_refined.status[i_ref,j_ref] = 1
        end
        if grid_refined.ref_level[i_ref, j_ref] < new_ref_level
            grid_refined.ref_level[i_ref, j_ref] = new_ref_level
        end
    end

    return n_calc
end
