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

#----------------------------------------------------------------------------------

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR/ANALYSYS",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_ABE0_BASE.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar",
    tparams = [TP("GAIN", 0.5), TP("NITS", 13), TP("GAIN_ONE_AFTER_NITS", 6)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = 32,
    nits = [13],
    gain = [0.5],
    tparams_local = [
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA",   0.0, flag=0), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=0), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=0), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=0), TP("OM2DOT", 0.0, flag=1)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=0)],
        [TP("P_DELTA", 100.0, flag=1), TP("IDOT", 0.0, flag=1), TP("I2DOT", 0.0, flag=1), TP("XOMDOT", 0.0, flag=1), TP("OM2DOT", 0.0, flag=1)]
        ],
    print_output = true
)

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

function extract_parameter(results_global_iters, param_name)
    N = length(results_global_iters.all_global_iterations)
    param_values = Vector{Float64}()
    param_errors = Vector{Float64}()
    for i in 1:N
        if results_global_iters.all_global_iterations[i].last_internal_iteration.error == TempoOutputError()
            try 
                if param_name == :chisqr
                    push!(param_values, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.chisqr[1])
                    push!(param_errors, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.chisqr[2])
                else
                    param_ind = results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters_order[param_name]
                    push!(param_values, results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].post_fit)
                    push!(param_errors, results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].uncertainty)
                end
            catch error
            
            end
        end
    end
    return param_values, param_errors
end