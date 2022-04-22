module GravityTools

using StructureSolver
using DelimitedFiles

include("Utils.jl")
include("MMDiagram.jl")
include("ClassicalTest.jl")

export read_DEFGrid, interpolate_DEFSmallGrid, interpolate_in_mass

end # module
