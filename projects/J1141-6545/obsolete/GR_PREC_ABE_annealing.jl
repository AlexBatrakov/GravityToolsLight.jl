using PyPlot
using Revise
using DelimitedFiles
using GravityToolsLight
using MultivariateStats
using Statistics
using GLM
using DataFrames
using JLD
using Printf
using Polynomials
using LsqFit
using Roots
using LsqFit
using Printf
using LombScargle
pygui(true)


using Distributed
addprocs(8)

@everywhere using Revise
@everywhere using GravityToolsLight

#----------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final/GR_PREC/ANNEALING",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PREC_ABE.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("NITS", 6)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

# results_basic = run_tempo_basic(basic_settings)

global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 3,
    nits = [4, 1, 4],
    gain = [1.0, 1.0, 1.0],
    tparams_local = [
        [TP("P_DELTA", flag=0), TP("P_PHI", flag=0), TP("C_DELTA", flag=0), TP("PREC_PHI", flag=0), TP("C_XI", flag=0)],
        [TP("P_DELTA", flag=1), TP("P_PHI", flag=1), TP("C_DELTA", flag=1), TP("PREC_PHI", flag=1), TP("C_XI", flag=1)],
        [TP("P_DELTA", flag=0), TP("P_PHI", flag=0), TP("C_DELTA", flag=0), TP("PREC_PHI", flag=0), TP("C_XI", flag=0)],
        ]
    )

# results = run_tempo_global_iters(basic_settings, global_iters_settings)

annealing_settings = AnnealingSettings(
    parameters = [AP("P_DELTA",  0.0, 0.0, 180.0, 0.0, true),
                  AP("P_PHI",    0.0, 0.0, 360.0, 0.0, true),
                  AP("C_DELTA",  0.0, 0.0, 180.0, 0.0, true),
                  AP("PREC_PHI", 0.0, 0.0, 360.0, 0.0, true),
                  AP("C_XI",     0.0, 0.0, 3.0,   0.0, false)
                ],
    energy_scale = 1000.0,
    quenching_factor = 5.0,
    initial_temperature = 1.0,
    cooling_rate = 0.95,
    minimum_temperature = 1e-5,
    step_law = 0.5,
    iterations_per_temp = 1,
    max_iterations = 100,
    parallel = true,
    gradient = true
)

# results_annealing = run_tempo_annealing(basic_settings, annealing_settings)

# @time results_annealing = run_tempo_annealing_parallel(basic_settings, annealing_settings)

@time results_annealing = run_tempo_annealing_parallel(basic_settings, global_iters_settings, annealing_settings)