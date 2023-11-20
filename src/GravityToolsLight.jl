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
export TestParameters
export G_CAV, M_sun, c, d, rad
# Mathematical constants and variables
export lvl_1σ, lvl_2σ, lvl_3σ, lvl_4σ, lvl_5σ, lvl_6σ, lvl_7σ, lvl_68CL, lvl_90CL, lvl_95CL, lvl_99CL
export LinRule, LogRule, RangeVariable, ValueVariable, Variable, Var
export get_label

include("AdaptiveGridFramework/SimpleGrid.jl")
include("AdaptiveGridFramework/Refinement2DGrid.jl")
export SimpleGrid, precalculate_Grid, refine_Grid, grid_size_counter

include("AdaptiveGridFramework/AdaptiveRefinement2DGrid.jl")

export FullUnit, DiffUnit, ContourUnit, DiffContourUnit, RefinementSettings, AdaptiveRefinement2DGrid, precalculate_2DGrid!, refine_2DGrid, calculate_2DGrid!

# GeneralTempoFramework
# Include and export from GeneralTempoFramework
# Включение и экспорт из Essentials
include("GeneralTempoFramework/Essentials/AbstractTempo.jl")
export AbstractTempoVersion, Tempo, Tempo2, get_tempo_directory, get_tempo2_directory, get_tempo_command
include("GeneralTempoFramework/Essentials/TempoParameters.jl")
export GeneralTempoParameter, TP, extract_GeneralTempoParameter, get_par_file_representation
include("GeneralTempoFramework/Essentials/TempoParFile.jl")
export TempoParFile, read_par_file!, write_par_file, update_par_file, update_par_file, modify_tparam!, extend_par_file!, generate_par_file_path
include("GeneralTempoFramework/Essentials/TempoOutputResults.jl")
export AbstractTempoResult, BasicTempoOutputResult, FitParameter, DetailedTempoOutputResult, TempoOutputError
export CalculatedResults, TempoRunErrorOutput, SingleTempoRunResult

# Включение и экспорт из Settings
include("GeneralTempoFramework/Settings/BasicTempoSettings.jl")
export BasicTempoSettings, BasicTempoKeys
include("GeneralTempoFramework/Settings/GlobalItersSettings.jl")
export GlobalIterationsSettings
include("GeneralTempoFramework/Settings/ParameterSweepSettings.jl")
export ParameterSweepSettings
include("GeneralTempoFramework/Settings/GeneralSettings.jl")
export GeneralSettings, GeneralTempoFramework


# Включение и экспорт из Execution
include("GeneralTempoFramework/Execution/RunBasic.jl")
export run_tempo_basic, format_and_validate_par_file
include("GeneralTempoFramework/Execution/RunGlobalIters.jl")
export run_tempo_global_iters
include("GeneralTempoFramework/Execution/RunParameterSweep.jl")
export run_tempo_parameter_sweep
include("GeneralTempoFramework/Execution/RunGeneral.jl")
export run_tempo_general

include("GeneralTempoFramework/Utils/TempoUtils.jl")
#export TempoUtils
include("GeneralTempoFramework/Utils/PlotUtils.jl")
export plot_fit_results, plot_fit_results_all
include("GeneralTempoFramework/Utils/EFACs_EQUADs.jl")
export calculate_EFACs_EQUADs


include("PhysicalFramework/PhysicalFramework.jl")
include("PostKeplerianFramework/PKFramework.jl")
export DEF, GR, Object, BinarySystem, Settings, DEFPhysicalFramework, read_grid!, interpolate_mgrid!,
    interpolate_psr!, interpolate_comp!, interpolate_bnsys!, calculate_PK_params!, calculate_X_params!
export read_DEFGrid, interpolate_DEFMassGrid, interpolate_NS
export PKFramework, find_initial_masses, check_terms_in_chisqr, find_best_masses, optimize_PK_method, find_masses, obs_params_dataset, ObsParams

include("Obsolete/TempoFramework.jl")
include("Obsolete/Tempo2Framework.jl")
include("Obsolete/EOSAgnosticFramework.jl")
export EOSAgnosticTest
export TempoFramework, GeneralTest, TempoSettings, GridSetttings, calculate!, cut_ddstg_grid!
export TempoParameter, get_TempoParameter, update_pf_theory!, modify_par_file, run_tempo, get_par_file_work, read_params
export calculate_t2!


Base.Float64(m::Measurement{Float64}) = Float64(m.val)


# include("Obsolete/ClassicalTest.jl")
# export ClassicalTest, find_initial_masses, find_best_masses, check_terms_in_chisqr, calculate!, ct_dataset
# include("Obsolete/MMDiagram.jl")
# export MMDiagram
# include("Obsolete/DDSTGTest.jl")
# export DDSTGTest, DDSTGTestSettings, cut_ddstg_grid!

end # module
