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
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_test",
    version = Tempo2(),
    par_file_init = "DDSTG_test.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("NITS", 1), TP("TNRedC", 84)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics   
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551

N_iters = 9

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = ones(N_iters),
    gain = ones(N_iters),
    tparams_local = [[TP("FINISH", time_finish)] for i in 1:N_iters],
    flags = ["" for i in 1:N_iters],
    fit_EFACs_EQUADs = [1 < i < N_iters - 2 ? true : false for i in 1:N_iters],
    print_output = true
)


results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)


#----------------------------------------------------------------------------------
basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_test",
    version = Tempo2(),
    par_file_init = "DDSTG_test.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("TNRedC", 84)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics   
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551
time_validation = time_finish - 200.0

N_iters = 9

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = ones(N_iters),
    gain = ones(N_iters),
    tparams_local = [
        if i < N_iters
            [TP("FINISH", time_validation, flag=1)]
        else
            [TP("FINISH", time_finish)]
        end
        for i in 1:N_iters], 
    flags = [i < N_iters ? "" : "-nofit" for i in 1:N_iters],
    fit_EFACs_EQUADs = [1 < i < N_iters - 2 ? true : false for i in 1:N_iters],
    print_output = true
)


results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)


#----------------------------------------------------------------------------------
basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_test",
    version = Tempo2(),
    par_file_init = "DDSTG_test.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("TNRedC", 30)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics   
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551
time_validation = time_finish - 200.0

N_iters = 13

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = ones(N_iters),
    gain = ones(N_iters),
    tparams_local = [
        if i in (1, 3, 5, 7, 9, 10, 11, 12)
            [TP("FINISH", time_validation, flag=1), TP("PL_T0", flag=0), TP("PL_PB", flag=0), TP("PL_ECC", flag=0), TP("PL_X", flag=0), TP("PL_OM", flag=0)]
        elseif i in (2, 4, 6, 8)
            [TP("FINISH", time_validation, flag=1), TP("PL_T0", flag=1), TP("PL_PB", flag=1), TP("PL_ECC", flag=1), TP("PL_X", flag=1), TP("PL_OM", flag=1)]
        else
            [TP("FINISH", time_finish)]
        end
        for i in 1:N_iters], 
    flags = [i < N_iters ? "" : "-nofit" for i in 1:N_iters],
    fit_EFACs_EQUADs = [i in (3, 5, 7, 9, 10) ? true : false for i in 1:N_iters],
    print_output = true
)


results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)