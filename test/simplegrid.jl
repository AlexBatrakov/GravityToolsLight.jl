using PyPlot
using GravityTools

function difference_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
    val1_cell = @view grid.value[:val1][i_cell:i_cell+1,j_cell:j_cell+1]
    val1_cell_min = minimum(val1_cell)
    val1_cell_max = maximum(val1_cell)
    return val1_cell_max - val1_cell_min > 0.1 ? true : false
end

function contour_cell_selector(i_cell::Int64, j_cell::Int64, grid::SimpleGrid)
    val1_cell = @view grid.value[:val1][i_cell:i_cell+1,j_cell:j_cell+1]
    val1_cell_min = minimum(val1_cell)
    val1_cell_max = maximum(val1_cell)
    return val1_cell_min < 0.5 < val1_cell_max
end

#-------------------------------------------------------------------------------------------------------------------------

x = [0.0, 1.0, 2.0]
y = [0.0, 1.0, 2.0, 3.0]
func(x,y) = sin(x^2*y)
func_full(x,y) = ((:val1,), (sin(x^2*y),))
dummy_params_func!(grid::SimpleGrid) = nothing

grid = SimpleGrid(Dict(), x, y)
precalculate_Grid(grid, func_full, dummy_params_func!)

grid_refined = refine_Grid(grid, func_full, contour_cell_selector, dummy_params_func!)

grid_array = Array{SimpleGrid}(undef, 10)
grid_array[1] = grid
for i in 2:10
    grid_array[i] = refine_Grid(grid_array[i-1], func_full, contour_cell_selector, calc_params!)
end

pygui(true)
imshow(grid_array[end].value[:val1])
imshow(grid_array[end].status)