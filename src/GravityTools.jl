module GravityTools

using StructureSolver
using DelimitedFiles
using Statistics
using NLsolve

include("SimpleGrid.jl")
include("Physics.jl")
include("Utils.jl")
include("ClassicalTest.jl")
include("MMDiagram.jl")

export SimpleGrid, precalculate_Grid, refine_Grid, grid_size_counter
export DEF, GR, Object, BinarySystem, Settings, DEFPhysicalFramework, read_grid!, interpolate_mgrid!,
    interpolate_psr!, interpolate_comp!, interpolate_bnsys!, calculate_PK_params!
export read_DEFGrid, interpolate_DEFMassGrid, interpolate_NS
export ClassicalTest, find_initial_masses

export G_CAV, M_sun, c, d, rad
end # module
