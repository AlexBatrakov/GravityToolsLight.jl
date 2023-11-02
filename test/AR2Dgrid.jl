#--------------------------------------------------------------------------------------------------------------
# Single_core routines
using Revise
using GravityToolsLight
using PyPlot


ref_sets = RefinementSettings(
    desired_refinement_level = 5,
    parallel = false,
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5,0.7])
    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

x = GravityToolsLight.Var(name = "var1", min = 0.0, max = 4.0, N = 5, range_rule=:lin)
y = GravityToolsLight.Var(name = "var2", min = 0.0, max = 3.0, N = 5, range_rule=:lin)
target_function(x,y) = (val1 = sin(x^2*y),)
params_function!(grid) = nothing

grid_init = AdaptiveRefinement2DGrid(x, y, ref_sets)

@time precalculate_2DGrid!(grid_init, target_function, params_function!)
grid_init.vars[:val1]

@time grid_refined = refine_2DGrid(grid_init, target_function, params_function!)
grid_refined.vars[:val1]

grid_array = Array{AdaptiveRefinement2DGrid}(undef, 9)
grid_array[1] = grid_init
for i in 2:9
    grid_array[i] = refine_2DGrid(grid_array[i-1], target_function, params_function!)
end

#--------------------------------------------------------------------------------------------------------------
# Parallel routines
using Revise
using GravityTests
using PyPlot
using Distributed

addprocs(8)

ref_sets = RefinementSettings(
    desired_refinement_level = 9,
    parallel = true,
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

x = GravityToolsLight.Var(name = "var1", min = 0.0, max = 4.0, N = 5, range_rule=:lin)
y = GravityToolsLight.Var(name = "var2", min = 0.0, max = 3.0, N = 5, range_rule=:lin)
@everywhere target_function(x,y) = (val1 = sin(x^2*y),)
@everywhere params_function!(grid) = nothing

grid_init = AdaptiveRefinement2DGrid(x, y, ref_sets)

@time precalculate_2DGrid!(grid_init, target_function, params_function!)
grid_init.vars[:val1]

@time grid_refined = refine_2DGrid(grid_init, target_function, params_function!)
grid_refined.vars[:val1]

grid_array = Array{AdaptiveRefinement2DGrid}(undef, 9)
grid_array[1] = grid_init
for i in 2:9
    grid_array[i] = refine_2DGrid(grid_array[i-1], target_function, params_function!)
end

grid_refined = calculate_2DGrid!(grid_init, target_function, params_function!)