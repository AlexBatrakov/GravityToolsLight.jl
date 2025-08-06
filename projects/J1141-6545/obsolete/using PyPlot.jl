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
pygui(true)


#----------------------------------------------------------------------------------
basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_RN/test",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PREC_ABE_BEND_RN.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("NITS", 1), TP("TNRedC", 83), TP("GAIN_ONE_AFTER_NITS", 5)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551
time_validation = time_finish - 200.0
# TNRedFlow = log10((time_validation - time_start) / (time_finish - time_start))


global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 8,
    nits = [1,10,10,10,10,10,10],
    gain = [1,0.2,0.2,0.2,0.2,0.2,0.2,0.2],
    tparams_local = [
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)],
        [TP("FINISH", time_finish)]
        ], 
    flags = ["", "", "", "", "", "", "", ""],
    fit_EFACs_EQUADs = [true, true, true, true, true, true, true, false],
    print_output = true
)


results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

