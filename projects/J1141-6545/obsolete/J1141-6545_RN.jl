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
using HypothesisTests
using KernelDensity
using QuadGK

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
    values = collect(-3e-14:5e-15:3e-14)
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


#-------------------------------------------------------------------------------------

xpbdot = vec(readdlm("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/original_data/J1141-6545_PbdotExt-Tot.dat"))
histogram = fit(Histogram, xpbdot, nbins=200)
pdf_values = histogram.weights / sum(histogram.weights)
bin_edges = histogram.edges[1]
plot(bin_edges[1:end-1], pdf_values)

pdf_mul = pdf_arr .* pdf_values 
pdf_mul ./= sum(pdf_mul)
plot(bin_edges[1:end-1], pdf_mul)

fig, ax = subplots()
hist(xpbdot, bins=200, density=true)

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

#-------------------------------------------------------------------------------------

xpbdot_prior = vec(readdlm("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/original_data/J1141-6545_PbdotExt-Tot.dat"))

xpbdot_prior .+= 10e-15

# Best fit 
chisqr_min, mu, sigma = 20914.19, -8.29481968690372e-15,  4.7427e-15 

xpbdot_normal_likelihood(xpbdot) = pdf(Normal(mu, sigma), xpbdot)
prior_kde = kde(xpbdot_prior)
xpbdot_prior_kde(xpbdot) = pdf(prior_kde, xpbdot)

delta_chisqr_normal(xpbdot) = -2 * log(xpbdot_normal_likelihood(xpbdot)) + 2 * log(xpbdot_normal_likelihood(mu))

xpbdot_fine_arr = collect(-5e-14:1e-16:5e-14)

fig, ax = subplots()
plot(xpbdot_fine_arr, xpbdot_normal_likelihood.(xpbdot_fine_arr), label = "Normal likelihood")
#hist(xpbdot, bins=200, density=true)
plot(xpbdot_fine_arr, xpbdot_prior_kde(xpbdot_fine_arr), label = "Prior")

product_kde_likelihood(xpbdot) = xpbdot_normal_likelihood(xpbdot) * xpbdot_prior_kde(xpbdot)
normalizing_constant_kde_likelihood, _ = quadgk(product_kde_likelihood, -5e-14, 5e-14)
normalized_posterior(xpbdot) = product_kde_likelihood(xpbdot) / normalizing_constant_kde_likelihood

plot(xpbdot_fine_arr, normalized_posterior.(xpbdot_fine_arr), label = "Posterior from Normal")
xlabel("XPBDOT, s/s")
ylabel("pdf")
legend()


objective_function(x) = -normalized_posterior(x)
result = optimize(objective_function, -5e-14, 5e-14)
max_posterior = -Optim.minimum(result)
max_posterior_value = Optim.minimizer(result)

# use sweep over XPBDOT to check likelihood 

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3), TP("XPBDOT", 5.933280313096283e-15, flag = 0)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )
bsets = basic_settings

results_basic = run_tempo_basic(basic_settings)
chisqr_min = 20914.19

global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 2,
    nits = [3,3],
    gain = [1.0, 1.0],
    tparams_local = [
        [TP("XPBDOT", -2e-14, flag = 0)],
        [TP("XPBDOT", -2.252291968690372e-14, flag = 0)],
        ]
    )
gisets = global_iters_settings

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)



parameter_sweep_settings = ParameterSweepSettings(
    parameter_name = "XPBDOT",
    #values = collect(-3e-14:1e-15:3e-14)
    values = collect(mu - 3*sigma : sigma : mu + 3.1*sigma)
)

results_sweep = run_tempo_parameter_sweep(basic_settings, parameter_sweep_settings)

chisqr_tempo_arr = [results_sweep[iter].last_internal_iteration.result.chisqr for iter in 1:length(parameter_sweep_settings.values)]

#--------------------------------------

xpbdot_tempo_range3 = mu - 1*sigma : sigma : mu + 1.1*sigma
xpbdot_tempo_arr3 = collect(xpbdot_tempo_range3)
chisqr_tempo_arr3 = [20915.51, 20914.19, 20914.87]

delta_chisqr_tempo_arr3 = chisqr_tempo_arr3 .- chisqr_min

#  pdf_tempo_arr5 = exp.(-0.5*delta_chisqr_tempo_arr5)

# bin_width5 = sigma

# pdf_tempo_arr_norm5 = pdf_tempo_arr5 ./ sum(pdf_tempo_arr5) ./ bin_width5

# plot(xpbdot_tempo_arr, pdf_tempo_arr_norm, ".", label = "true tempo pdf")
# legend()

delta_chisqr_interp_cubic3 = cubic_spline_interpolation(xpbdot_tempo_range3, delta_chisqr_tempo_arr3, extrapolation_bc = Line())

objective_function(x) = delta_chisqr_interp_cubic3(x)
result = optimize(objective_function, -5e-14, 5e-14)
delta_chisqr_interp_cubic3_min = Optim.minimum(result)
max_posterior_value = Optim.minimizer(result)

# Min value is not important -- we anyway normilize the distribution

likelihood_interp_cubic3(xpbdot) = exp(-0.5 * delta_chisqr_interp_cubic3(xpbdot))
normalizing_constant_likelihood3, _ = quadgk(likelihood_interp_cubic3, -5e-14, 5e-14)

plot(xpbdot_tempo_arr3, likelihood_interp_cubic3.(xpbdot_tempo_arr3) ./ normalizing_constant_likelihood3, "d", color = "grey")
plot(xpbdot_fine_arr, likelihood_interp_cubic3.(xpbdot_fine_arr) ./ normalizing_constant_likelihood3, label = "Splines likelihood 3", color = "grey")
xlabel("XPBDOT, s/s")
ylabel("pdf")
legend()

#--------------------------------------

xpbdot_tempo_range5 = mu - 2*sigma : sigma : mu + 2.1*sigma
xpbdot_tempo_arr5 = collect(xpbdot_tempo_range5)
chisqr_tempo_arr5 = [20919.48, 20915.51, 20914.19, 20914.87, 20916.72]

delta_chisqr_tempo_arr5 = chisqr_tempo_arr5 .- chisqr_min

#  pdf_tempo_arr5 = exp.(-0.5*delta_chisqr_tempo_arr5)

# bin_width5 = sigma

# pdf_tempo_arr_norm5 = pdf_tempo_arr5 ./ sum(pdf_tempo_arr5) ./ bin_width5

# plot(xpbdot_tempo_arr, pdf_tempo_arr_norm, ".", label = "true tempo pdf")
# legend()

delta_chisqr_interp_cubic5 = cubic_spline_interpolation(xpbdot_tempo_range5, delta_chisqr_tempo_arr5, extrapolation_bc = Line())

objective_function(x) = delta_chisqr_interp_cubic5(x)
result = optimize(objective_function, -5e-14, 5e-14)
delta_chisqr_interp_cubic5_min = Optim.minimum(result)
max_posterior_value = Optim.minimizer(result)

likelihood_interp_cubic5(xpbdot) = exp(-0.5 * delta_chisqr_interp_cubic5(xpbdot))
normalizing_constant_likelihood5, _ = quadgk(likelihood_interp_cubic5, -5e-14, 5e-14)

plot(xpbdot_tempo_arr5, likelihood_interp_cubic5.(xpbdot_tempo_arr5) ./ normalizing_constant_likelihood5, "d", color = "violet")
plot(xpbdot_fine_arr, likelihood_interp_cubic5.(xpbdot_fine_arr) ./ normalizing_constant_likelihood5, label = "Splines likelihood 5", color = "violet")
legend()

#--------------------------------------

xpbdot_tempo_range7 = mu - 3*sigma : sigma : mu + 3.1*sigma
xpbdot_tempo_arr7 = collect(xpbdot_tempo_range7)
chisqr_tempo_arr7 = [20926.38, 20919.48, 20915.51, 20914.19, 20914.87, 20916.72, 20918.39]

delta_chisqr_tempo_arr7 = chisqr_tempo_arr7 .- chisqr_min

delta_chisqr_interp_cubic7 = cubic_spline_interpolation(xpbdot_tempo_range7, delta_chisqr_tempo_arr7, extrapolation_bc = Line())

objective_function(x) = delta_chisqr_interp_cubic7(x)
result = optimize(objective_function, -5e-14, 5e-14)
delta_chisqr_interp_cubic7_min = Optim.minimum(result)
max_posterior_value = Optim.minimizer(result)

likelihood_interp_cubic7(xpbdot) = exp(-0.5 * delta_chisqr_interp_cubic7(xpbdot))
normalizing_constant_likelihood7, _ = quadgk(likelihood_interp_cubic7, -5e-14, 5e-14)

plot(xpbdot_tempo_arr7, likelihood_interp_cubic7.(xpbdot_tempo_arr7) ./ normalizing_constant_likelihood7, "d", color = "cyan")
plot(xpbdot_fine_arr, likelihood_interp_cubic5.(xpbdot_fine_arr) ./ normalizing_constant_likelihood7, label = "Splines likelihood 7", color = "cyan")
legend()

#--------------------------------------

xpbdot_tempo_range51 = -2.1e-14:1e-15:3e-14
xpbdot_tempo_arr51 = collect(xpbdot_tempo_range51)
chisqr_tempo_arr51 = [20924.39, 20921.83, 20920.92, 20919.73, 20918.63, 20917.65, 20916.8, 20916.08, 20915.49, 20915.01, 20914.65, 20914.39, 20914.24, 20914.18, 20914.21, 20914.33, 20914.51, 20914.75, 20915.03, 20915.36, 20915.71, 20916.08, 20916.46, 20916.85, 20917.24, 20917.64, 20918.03, 20918.42, 20918.8, 20919.18, 20919.55, 20919.93, 20920.29, 20920.65, 20921.01, 20921.36, 20921.71, 20922.06, 20922.4, 20922.75, 20923.09, 20923.42, 20923.76, 20924.1, 20924.44, 20924.77, 20925.11, 20925.45, 20925.79, 20926.14, 20926.48]

delta_chisqr_tempo_arr51 = chisqr_tempo_arr51 .- chisqr_min

delta_chisqr_interp_cubic51 = cubic_spline_interpolation(xpbdot_tempo_range51, delta_chisqr_tempo_arr51, extrapolation_bc = Line())

objective_function(x) = delta_chisqr_interp_cubic51(x)
result = optimize(objective_function, -5e-14, 5e-14)
delta_chisqr_interp_cubic51_min = Optim.minimum(result)
xpbdot_cubic51 = Optim.minimizer(result)

likelihood_interp_cubic51(xpbdot) = exp(-0.5 * delta_chisqr_interp_cubic51(xpbdot))
normalizing_constant_likelihood51, _ = quadgk(likelihood_interp_cubic51, -5e-14, 5e-14)

plot(xpbdot_tempo_arr51, likelihood_interp_cubic51.(xpbdot_tempo_arr51) ./ normalizing_constant_likelihood51, "d", color = "red")
plot(xpbdot_fine_arr, likelihood_interp_cubic51.(xpbdot_fine_arr) ./ normalizing_constant_likelihood51, label = "Splines likelihood 51", color = "red")
legend()

#--------------------------------------

fig, ax = subplots()


plot(xpbdot_fine_arr, delta_chisqr_interp_cubic51.(xpbdot_fine_arr) .- delta_chisqr_interp_cubic51_min, label = "Splines likelihood 51", color = "c")


plot(xpbdot_fine_arr, delta_chisqr_interp_cubic7.(xpbdot_fine_arr) .- delta_chisqr_interp_cubic7_min, label = "Splines likelihood 7", color = "r")


plot(xpbdot_fine_arr, delta_chisqr_interp_cubic5.(xpbdot_fine_arr) .- delta_chisqr_interp_cubic5_min, label = "Splines likelihood 5", color = "g")


plot(xpbdot_fine_arr, delta_chisqr_interp_cubic3.(xpbdot_fine_arr) .- delta_chisqr_interp_cubic3_min, label = "Splines likelihood 3", color = "b")

plot(xpbdot_tempo_arr51, delta_chisqr_interp_cubic51.(xpbdot_tempo_arr51) .- delta_chisqr_interp_cubic51_min, "d", color = "c")
plot(xpbdot_tempo_arr7, delta_chisqr_interp_cubic7.(xpbdot_tempo_arr7) .- delta_chisqr_interp_cubic7_min, "d", color = "r")
plot(xpbdot_tempo_arr5, delta_chisqr_interp_cubic5.(xpbdot_tempo_arr5) .- delta_chisqr_interp_cubic5_min, "d", color = "g")
plot(xpbdot_tempo_arr3, delta_chisqr_interp_cubic3.(xpbdot_tempo_arr3) .- delta_chisqr_interp_cubic3_min, "d", color = "b")

legend()

xpbdot_prior_chisqr(xpbdot) = -2 * log(xpbdot_prior_kde(xpbdot))
plot(xpbdot_fine_arr, xpbdot_prior_chisqr.(xpbdot_fine_arr) .- minimum(xpbdot_prior_chisqr.(xpbdot_fine_arr)), label = "Prior")

xpbdot_normal_likelihood_chisqr(xpbdot) = -2 *  log(xpbdot_normal_likelihood(xpbdot))
plot(xpbdot_fine_arr, xpbdot_normal_likelihood_chisqr.(xpbdot_fine_arr) .- minimum(xpbdot_normal_likelihood_chisqr.(xpbdot_fine_arr)), label = "Normal")


xlabel("XPBDOT, s/s")
ylabel(L"\Delta \chi^2")

legend()

#--------------------------------------

function get_chisqr_posterior(xpbdot, prior::Function, likelihood::Function)
    product_prior_likelihood(xpbdot) = prior(xpbdot) * likelihood(xpbdot)
    normalizing_constant, _ = quadgk(product_prior_likelihood, -5e-14, 5e-14)
    posterior(xpbdot) = product_prior_likelihood(xpbdot) / normalizing_constant

    result = optimize(x -> -prior(x), -5e-14, 5e-14)
    max_prior_prob = -Optim.minimum(result)
    max_prior_xpbdot = Optim.minimizer(result)

    result = optimize(x -> -posterior(x), -5e-14, 5e-14)
    max_posterior_prob = -Optim.minimum(result)
    max_posterior_xpbdot = Optim.minimizer(result)

    delta_chisqr_posterior(xpbdot) = -2 * log(posterior(xpbdot)) + 2 * log(posterior(max_prior_xpbdot))

    return delta_chisqr_posterior(xpbdot)
end


function get_delta_chisqr_posterior(xpbdot_arr, prior::Function, posterior::Function)
    result = optimize(x -> -prior(x), -5e-14, 5e-14)
    max_prior_prob = -Optim.minimum(result)
    max_prior_xpbdot = Optim.minimizer(result)

    result = optimize(x -> -posterior(x), -5e-14, 5e-14)
    max_posterior_prob = -Optim.minimum(result)
    max_posterior_xpbdot = Optim.minimizer(result)

    delta_chisqr_posterior(xpbdot) = -2 * log(posterior(xpbdot)) + 2 * log(posterior(max_prior_xpbdot))

    return delta_chisqr_posterior.(xpbdot_arr), max_posterior_xpbdot, delta_chisqr_posterior(max_posterior_xpbdot)
end

fig, ax = subplots()

product_kde_normal_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * xpbdot_normal_likelihood(xpbdot)
normalizing_constant_kde_normal_likelihood, _ = quadgk(product_kde_normal_likelihood, -5e-14, 5e-14)
normal_posterior(xpbdot) = product_kde_normal_likelihood(xpbdot) / normalizing_constant_kde_normal_likelihood
plot(xpbdot_fine_arr, normal_posterior.(xpbdot_fine_arr), label = "Normal posterior")


product_kde_cubic3_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_interp_cubic3(xpbdot)
normalizing_constant_kde_cubic3_likelihood, _ = quadgk(product_kde_cubic3_likelihood, -5e-14, 5e-14)
cubic3_posterior(xpbdot) = product_kde_cubic3_likelihood(xpbdot) / normalizing_constant_kde_cubic3_likelihood
plot(xpbdot_fine_arr, cubic3_posterior.(xpbdot_fine_arr), label = "Cubic3 posterior")

product_kde_cubic5_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_interp_cubic5(xpbdot)
normalizing_constant_kde_cubic5_likelihood, _ = quadgk(product_kde_cubic5_likelihood, -5e-14, 5e-14)
cubic5_posterior(xpbdot) = product_kde_cubic5_likelihood(xpbdot) / normalizing_constant_kde_cubic5_likelihood
plot(xpbdot_fine_arr, cubic5_posterior.(xpbdot_fine_arr), label = "Cubic5 posterior")

product_kde_cubic7_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_interp_cubic7(xpbdot)
normalizing_constant_kde_cubic7_likelihood, _ = quadgk(product_kde_cubic7_likelihood, -5e-14, 5e-14)
cubic7_posterior(xpbdot) = product_kde_cubic7_likelihood(xpbdot) / normalizing_constant_kde_cubic7_likelihood
plot(xpbdot_fine_arr, cubic7_posterior.(xpbdot_fine_arr), label = "Cubic7 posterior")

product_kde_cubic51_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_interp_cubic51(xpbdot)
normalizing_constant_kde_cubic51_likelihood, _ = quadgk(product_kde_cubic51_likelihood, -5e-14, 5e-14)
cubic51_posterior(xpbdot) = product_kde_cubic51_likelihood(xpbdot) / normalizing_constant_kde_cubic51_likelihood
plot(xpbdot_fine_arr, cubic51_posterior.(xpbdot_fine_arr), label = "Cubic51 posterior")

xlabel("XPBDOT, s/s")
ylabel("pdf")

legend()

#--------------------------------------

prior = xpbdot_prior_kde
posteriors = [normal_posterior, cubic3_posterior, cubic5_posterior, cubic7_posterior, cubic51_posterior]
delta_chisqr_tempo_arr = [delta_chisqr_normal, delta_chisqr_interp_cubic3, delta_chisqr_interp_cubic5, delta_chisqr_interp_cubic7, delta_chisqr_interp_cubic51]
labels = ["Normal", "Cubic3", "Cubic5", "Cubic7", "Cubic51"]

fig, ax = subplots()
for i in eachindex(posteriors)
    delta_chisqr_tempo = delta_chisqr_tempo_arr[i]


    posterior = posteriors[i]
    delta_chisqr_posterior_fine_arr, max_posterior_xpbdot, posterior_correction = get_delta_chisqr_posterior(xpbdot_fine_arr, prior, posterior)
    plot(xpbdot_fine_arr, delta_chisqr_tempo.(xpbdot_fine_arr) + delta_chisqr_posterior_fine_arr, label = labels[i] * " " * string(delta_chisqr_tempo(max_posterior_xpbdot) + posterior_correction))
    plot(max_posterior_xpbdot, delta_chisqr_tempo(max_posterior_xpbdot) + posterior_correction, "d", color="red")
    println(labels[i], " ", delta_chisqr_tempo(max_posterior_xpbdot) + posterior_correction)
end
xlabel("XPBDOT, s/s")
ylabel(L"\Delta \chi^2 \mathrm{, \, posterior}")
legend()

#--------------------------------------
likelihoods = [xpbdot_normal_likelihood, likelihood_interp_cubic3, likelihood_interp_cubic5, likelihood_interp_cubic7, likelihood_interp_cubic51]
labels = ["Normal", "Cubic3", "Cubic5", "Cubic7", "Cubic51"]

for i in eachindex(likelihoods)
    likelihood = likelihoods[i]

    prior = xpbdot_prior_kde

    result = optimize(x -> - likelihood(x), -5e-14, 5e-14)
    xpbdot_cubic51 = Optim.minimizer(result)
    normalizing_constant_likelihood, _ = quadgk(likelihood, -5e-14, 5e-14)

    P, _ = quadgk(x -> prior(x) * likelihood(x) / normalizing_constant_likelihood, -5e-14, 5e-14)

    P_delta = likelihood(xpbdot_cubic51) / normalizing_constant_likelihood

    delta_chisqr_margin = -2 * log(P) + 2 * log(P_delta)
    println(labels[i], " ", delta_chisqr_margin)
end

#-------------------------------------------------------------------------------------



chisqr_tempo_arr = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 20924.39, 20921.83, 20920.92, 20919.73, 20918.63, 20917.65, 20916.8, 20916.08, 20915.49, 20915.01, 20914.65, 20914.39, 20914.24, 20914.18, 20914.21, 20914.33, 20914.51, 20914.75, 20915.03, 20915.36, 20915.71, 20916.08, 20916.46, 20916.85, 20917.24, 20917.64, 20918.03, 20918.42, 20918.8, 20919.18, 20919.55, 20919.93, 20920.29, 20920.65, 20921.01, 20921.36, 20921.71, 20922.06, 20922.4, 20922.75, 20923.09, 20923.42, 20923.76, 20924.1, 20924.44, 20924.77, 20925.11, 20925.45, 20925.79, 20926.14, 20926.48]
#chisqr_tempo_arr = map(chisqr -> isnan(chisqr) ? Inf : chisqr, chisqr_tempo_arr)
chisqr_tempo_arr = chisqr_tempo_arr[10:end]


chisqr_tempo_min = minimum(chisqr_tempo_arr)
delta_chisqr_tempo_arr = chisqr_tempo_arr .- chisqr_tempo_min

bin_width = 1e-15
#pdf_chisqr_values = pdf.(Chisq(1), delta_chisqr)

pdf_tempo_arr = exp.(-0.5*delta_chisqr_tempo_arr)
#pdf_chisqr_arr = pdf.(Chisq(1), delta_chisqr_arr)
xpbdot_tempo_arr = collect(-3e-14:1e-15:3e-14)[10:end]
xpbdot_tempo_arr = collect(mu - 2*sigma : sigma : mu + 2.1*sigma)

#kde_result = kde(xpbdot_chisqr_arr, weights=pdf_chisqr_arr)

#plot(kde_result.x, kde_result.density, label="KDE of xpbdot")

pdf_tempo_arr_norm = pdf_tempo_arr ./ sum(pdf_tempo_arr) ./ bin_width

plot(xpbdot_tempo_arr, pdf_tempo_arr_norm, label = "true tempo pdf")
legend()


log_likelihood_max = maximum(2 .* log.(xpbdot_normal_likelihood.(xpbdot_fine_arr)))

fig, ax = subplots()
plot(xpbdot_chisqr_arr, delta_chisqr_tempo_arr, ".", label = "tempo output")
plot(mu, chisqr_min - chisqr_tempo_min, "x", color="red")
#ax.errorbar(collect(-3e-14:1e-15:3e-14), chisqr_tempo_arr, yerr=0.005)
plot(xpbdot_fine_arr, log_likelihood_max .- 2 .* log.(xpbdot_normal_likelihood.(xpbdot_fine_arr)), label = "normal likelihood")


delta_chisqr_interp_cubic = cubic_spline_interpolation(-2.1e-14:1e-15:3e-14, delta_chisqr_tempo_arr, extrapolation_bc = Line())
delta_chisqr_interp_cubic = cubic_spline_interpolation(mu - 2*sigma : sigma : mu + 2.1*sigma, delta_chisqr_tempo_arr, extrapolation_bc = Line())

pdf_iterpol(xpbdot) = exp(-0.5 * delta_chisqr_interp_cubic(xpbdot))

normalizing_constant, _ = quadgk(pdf_iterpol, -5e-14, 5e-14)

plot(xpbdot_fine_arr, delta_chisqr_interp_cubic.(xpbdot_fine_arr), label = "cubic cpline 5 points")



plot(xpbdot_fine_arr, pdf_iterpol.(xpbdot_fine_arr) .* normalizing_constant, label = "cubic cplines")


fitted_model = fit_mle(Normal, xpbdot_chisqr_arr, pdf_chisqr_arr)
dist = Normal(fitted_model.μ, fitted_model.σ)
plot(xpbdot_fine_arr, pdf.(dist, xpbdot_fine_arr), label = "inter norm")
legend()

xpbdot_values, pdf_tempo_arr = ([-3.0e-14, -2.4999999999999998e-14, -1.9999999999999997e-14, -1.5e-14, -9.999999999999998e-15, -4.999999999999998e-15, 0.0, 5.000000000000004e-15, 1.0000000000000002e-14, 1.5e-14, 2.0000000000000003e-14, 2.5e-14, 3.0e-14], [0.0, 0.0, 1.6356057394250625e12, 2.0226896296365285e13, 6.7492279838020125e13, 6.3561835442222984e13, 2.899180542680611e13, 1.0935487322998014e13, 4.2292015592668984e12, 1.7367459498631523e12, 7.386088596573391e11, 3.172748705050838e11, 1.3425869486994722e11])


pdf_SEP(x, μ, σ) = pdf(Normal(μ,σ), x)

function loss_function(params; x = xpbdot_chisqr_arr, observed_pdf = pdf_chisqr_arr)
    μ, σ = params
    predicted_pdf = pdf_SEP.(x, μ, σ)
    sum((predicted_pdf - observed_pdf).^2)
end


lower = [-Inf, 0]
upper = [Inf, Inf]
initial_params = [mu, sigma]
inner_optimizer = GradientDescent(linesearch=LineSearches.BackTracking(order=3))
result = optimize(loss_function, lower, upper, initial_params, Fminbox(inner_optimizer))
optimal_params = Optim.minimizer(result)




plot(xpbdot_chisqr_arr, pdf_chisqr_arr, label = "true tempo")
plot(xpbdot_fine_arr, pdf_SEP.(xpbdot_fine_arr, optimal_params...), label = "interpol")
legend()