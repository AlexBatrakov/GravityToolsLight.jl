module GravityTools

using StructureSolver
using DelimitedFiles

include("Physics.jl")
include("Utils.jl")
include("MMDiagram.jl")
include("ClassicalTest.jl")

export DEF, GR, Object, BinarySystem, Settings, DEFPhysicalFramework, reinitialize!
export read_DEFGrid, interpolate_DEFMassGrid, interpolate_NS

end # module
