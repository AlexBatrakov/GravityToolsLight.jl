using PyPlot
using Revise
using GravityToolsLight
using MultivariateStats
using Statistics
using GLM
using DataFrames
using JLD
pygui(true)

using Distributed
addprocs(8)

@everywhere using GravityToolsLight

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_fake.par",
    tim_file = "DDSTG_GR_fake.simulate",
    # par_file_init = "DDSTG_GR_ABE_test.par",
    # tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

# results_basic = run_tempo_basic(basic_settings)

annealing_settings = AnnealingSettings(
    parameters = [AP("P_DELTA", 0.0, 0.0, 180.0, 180.0, true), AP("P_PHI", 0.0, 0.0, 360.0, 360.0, true)],
    initial_temperature = 6.0,
    cooling_rate = 0.95,
    minimum_temperature = 0.0,
    iterations_per_temp = 1,
    max_iterations = 5,
    parallel = true
)

# results_annealing = run_tempo_annealing(basic_settings, annealing_settings)

@time results_annealing = run_tempo_annealing_parallel(basic_settings, annealing_settings)