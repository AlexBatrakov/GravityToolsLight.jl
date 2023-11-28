using Contour
using JLD
using ColorSchemes
using Statistics
using DelimitedFiles
using Distributions
using StructArrays
using StatsBase
using Interpolations
using QuadGK
using Optim
using Roots

using Revise
using GravityToolsLight
using PyPlot
pygui(true)

# using Distributed
# addprocs(8)
# @everywhere using GravityToolsLight
# #pkill -f /opt/julia-1.8.2/bin/julia

#-------------------------------------------------------------------------------------

# tsets = GeneralTempoSettings(
#     work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/model_WN",
#     version = Tempo2(),
#     par_file_init = "J1141-6545_full_DDSTG_RN.par",
#     tim_file = "J1141-6545_full.tim",
#     flags = "-nobs 34000 -newpar -writeres -residuals",
#     tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
#     keys = TempoKeys(silent=true, print_output=true, iterative_mode=true, fit_EFACs_EQUADs=true),
#     iters = 3,
#     nits = [3,3,3],
#     tparams_local = [
#         [TP("M2", 1.0, flag=0)],
#         [TP("M2", flag=1)],
#         [TP("M2", flag=1)]
#         ]
#     )

# parsed_results, output, stderr_output = run_tempo_single(tsets)

# GravityToolsLight.plot_iterations_data(parsed_results)

# plot_fit_results(parsed_results)

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/model_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_full_DDSTG_RN.par",
    tim_file = "J1141-6545_full.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )
bsets = basic_settings


results_basic = run_tempo_basic(basic_settings)
lii = results_basic.last_internal_iteration
result = results_basic.last_internal_iteration.result

result.fit_parameters 
result.chisqr
result.pre_post
result.XPBDOT

aii = results_basic.all_internal_iterations

extract_internal_iterations_values(results_basic, :MTOT, :post_fit)
extract_internal_iterations_values(results_basic, :chisqr, :post_fit)

#plot_fit_results(parsed_output)

global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 3,
    nits = [2,2,2],
    gain = [1.0, 1.0, 0.1],
    tparams_local = [
        [TP("MTOT", flag=0), TP("M2", 1.0, flag=0)],
        [TP("MTOT", flag=1), TP("M2", 1.0, flag=0)],
        [TP("MTOT", flag=1), TP("M2", flag=1)]
        ]
    )
gisets = global_iters_settings

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

extract_internal_iterations_values(results_global_iters, :MTOT, :post_fit)
extract_internal_iterations_values(results_global_iters, :M2, :post_fit)
extract_internal_iterations_values(results_global_iters, :chisqr)
plot_fit_results(parsed_output)

parameter_sweep_settings = ParameterSweepSettings(
    parameter_name = "XPBDOT",
    values = [-2e-14, -1e-14, 0.0, 1e-14, 2e-14]
)

sweep_results = run_tempo_parameter_sweep(basic_settings, global_iters_settings, parameter_sweep_settings)


general_settings = NewGeneralTempoSettings(basic_settings, global_iter_settings, parameter_sweep_settings)

# test = GeneralTest(
#     psrname = "J1141-6545",
#     eosname = "BSk22",
#     param1 = (name = "TNRedAmp", min = 8.0, max = 12.0, N = 5),
#     param2 = (name = "TNRedGam", min = 1.0, max = 4.0, N = 4)
#     )

test_params = TestParameters(
    Var(name = "TNRedAmp", min = 8.0, max = 12.0, N = 5, range_rule=:lin),
    Var(name = "TNRedGam", min = 1.0, max = 4.0, N = 4, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

# gsets = GridSetttings(
#     N_refinement = 0,
#     CL = [0.90],
#     refinement_type = "full",
#     delta_chisqr_max = 10,
#     delta_chisqr_diff = 1.0,
#     gr_in_chisqr = true
#     )

ref_sets = RefinementSettings(
    desired_refinement_level = 2,
    parallel = false,
    FullUnit(:chisqr)
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
#    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )



tf = GeneralTempoFramework(tsets, test_params, ref_sets)

calculate!(tf)



#-------------------------------------------------------------------------------------

par_file = TempoParFile(tf.tsets.par_file_init)
par_file.tparams[:TNRedAmp].value = -8.0
par_file.tparams[:TNRedGam].value = 3.0
tparam = GeneralTempoParameter("TNEF -be CASPSR", 0.991577)
extend_par_file!(par_file, tparam)
write_par_file(par_file, "test.par")

#-------------------------------------------------------------------------------------


basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_full_DDSTG_best2_WN.par",
#     par_file_init = "J1141-6545_until_2018_DDSTG_best.par",
    tim_file = "J1141-6545_full_WN.tim",
#     tim_file  = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )
bsets = basic_settings

# results_basic = run_tempo_basic(basic_settings)

parameter_sweep_settings = ParameterSweepSettings(
    parameter_name = "XPBDOT",
#    values = collect(-1.5e-14:1e-15:1.5e-14)
    values = collect(-1.16e-14:2.0e-16:1.82e-14)
    #values = [-1e-14, 0, 1e-14]
)

results_sweep = run_tempo_parameter_sweep(basic_settings, parameter_sweep_settings)

save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN/results_sweep.jld",  "results_sweep", results_sweep)
results_sweep = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN/results_sweep.jld", "results_sweep")

chisqr_arr = [results_sweep[iter].last_internal_iteration.result.chisqr for iter in 1:length(parameter_sweep_settings.values)]
pdf_arr = exp.(-0.5*(chisqr_arr .- minimum(filter(x -> !isnan(x), chisqr_arr))))
pdf_arr ./= sum(pdf_arr)
plot(parameter_sweep_settings.values, pdf_arr)


xpbdot = vec(readdlm("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/original_data/J1141-6545_PbdotExt-Tot.dat"))
histogram = fit(Histogram, xpbdot, nbins=200)
pdf_values = histogram.weights / sum(histogram.weights)
bin_edges = histogram.edges[1]
plot(bin_edges[1:end-1], pdf_values)

pdf_mul = pdf_arr .* pdf_values 
pdf_mul ./= sum(pdf_mul)
plot(bin_edges[1:end-1], pdf_mul)

#axvline(x=-4.33365376107607e-15)
#axvline(x=-4.33365376107607e-15+2.6574e-15, color="red")
#axvline(x=-4.33365376107607e-15-2.6574e-15, color="red")

xpbdot_ticks = collect(bin_edges[1:end-1])
pdf_mul_function = interpolate((xpbdot_ticks,), pdf_mul, Gridded(Linear()))
pdf_prior_function = interpolate((xpbdot_ticks,), pdf_values, Gridded(Linear()))

cdf = cumsum(pdf_mul) / sum(pdf_mul)

# Разделение CDF на равные части
n_points = 6
targets = [(i/n_points) * cdf[end] for i in 1:n_points-1]

# Вычисление центров сегментов
selected_points = []
for i in 1:length(targets)-1
    left_idx = findfirst(x -> x >= targets[i], cdf)
    right_idx = findfirst(x -> x >= targets[i+1], cdf)
    center_point = (xpbdot_ticks[left_idx] + xpbdot_ticks[right_idx]) / 2
    push!(selected_points, center_point)
end

selected_points # содержит выбранные точки для интегрирования

plot(selected_points, pdf_function.(selected_points), "o", color="violet")



integral_mul, error = quadgk(pdf_mul_function, xpbdot_ticks[1], xpbdot_ticks[end])
integral_prior, error = quadgk(pdf_prior_function, xpbdot_ticks[1], xpbdot_ticks[end])

integral_mul / integral_prior


XDOT_arr = [results_sweep[iter].last_internal_iteration.result.XDOT.post_fit for iter in 1:length(parameter_sweep_settings.values)]

plot(xpbdot_ticks, XDOT_arr)