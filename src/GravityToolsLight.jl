module GravityToolsLight

#using StructureSolver
using DelimitedFiles
using Distributions
using NLsolve
using Measurements
using Optim
using StructArrays
using Printf
using Distributed
using ColorSchemes
using PyPlot
using ProgressMeter

include("Utils.jl")
include("SimpleGrid.jl")

include("Refinement2DGrid.jl")

include("PhysicalFramework.jl")
#include("ClassicalTest.jl")
#include("MMDiagram.jl")
#include("DDSTGTest.jl")
include("PKFramework.jl")
include("TempoFramework.jl")
include("Tempo2Framework.jl")
include("EOSAgnosticFramework.jl")

Base.Float64(m::Measurement{Float64}) = Float64(m.val)

export get_label
export SimpleGrid, precalculate_Grid, refine_Grid, grid_size_counter
export DEF, GR, Object, BinarySystem, Settings, DEFPhysicalFramework, read_grid!, interpolate_mgrid!,
    interpolate_psr!, interpolate_comp!, interpolate_bnsys!, calculate_PK_params!, calculate_X_params!
export read_DEFGrid, interpolate_DEFMassGrid, interpolate_NS
#export ClassicalTest, find_initial_masses, find_best_masses, check_terms_in_chisqr, calculate!, ct_dataset
#export MMDiagram
#export DDSTGTest, DDSTGTestSettings, cut_ddstg_grid!

export PKFramework, find_initial_masses, check_terms_in_chisqr, find_best_masses, optimize_PK_method, find_masses, obs_params_dataset, ObsParams
export TempoFramework, GeneralTest, TempoSettings, GridSetttings, calculate!, cut_ddstg_grid!
export TempoParameter, get_TempoParameter, update_pf_theory!, modify_par_file, run_tempo, get_par_file_work, read_params
export calculate_t2!
export EOSAgnosticTest

export G_CAV, M_sun, c, d, rad
export  lvl_1σ, lvl_2σ, lvl_3σ, lvl_4σ, lvl_5σ, lvl_6σ, lvl_7σ, lvl_68CL, lvl_90CL, lvl_95CL, lvl_99CL
end # module
