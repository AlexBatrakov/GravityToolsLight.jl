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

#------------------------------------------------------------------------------------- 
# WE HAVE PRIOR INFORMATION ABOUT XPBDOT
# use sweep or global iterations over XPBDOT to check likelihood
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# FIRSTLY WE FIT FOR XPBDOT IN TEMPO
#-------------------------------------------------------------------------------------
basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

results_basic = run_tempo_basic(basic_settings)

# BEST FIT 
chisqr_min = results_basic.last_internal_iteration.result.chisqr
xpbdot_mu = results_basic.last_internal_iteration.result.XPBDOT.post_fit
xpbdot_sigma = results_basic.last_internal_iteration.result.XPBDOT.uncertainty

# WE OBTAINED MU AND SIGMA FOR NORMAL LIKELIHOOD
chisqr_min, xpbdot_mu, xpbdot_sigma = 20914.19, -8.29481968690372e-15,  4.7427e-15
# BUILD NORMAL LIKELIHOOD
xpbdot_normal_likelihood(xpbdot) = pdf(Normal(xpbdot_mu, xpbdot_sigma), xpbdot)
# BUILD NORMAL LIKELIHOOD IN CHISQR REPRESENTATION
delta_chisqr_normal_likelihood(xpbdot) = -2 * log(xpbdot_normal_likelihood(xpbdot)) + 2 * log(xpbdot_normal_likelihood(xpbdot_mu))


#-------------------------------------------------------------------------------------
# NOW WE UPLOAD PRIOR INFORMATION
#-------------------------------------------------------------------------------------
xpbdot_prior = vec(readdlm("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/original_data/J1141-6545_PbdotExt-Tot.dat"))

# IF WE WANT TO SHIFT PRIOR FOR INVESTIGATION REASONS
xpbdot_prior .+= 1e-14

# WE OBTAINED MU AND SIGMA FROM RAW PRIOR INFORMATION
xpbdot_prior_mu, xpbdot_prior_sigma = mean(xpbdot_prior), std(xpbdot_prior)

# BUILD TRUE PRIOR USING KDE
prior_kde = kde(xpbdot_prior)
xpbdot_prior_kde(xpbdot) = pdf(prior_kde, xpbdot)

# FIND MAXIMIMUM OF PRIOR
result = optimize(x -> -xpbdot_prior_kde(x), -5e-14, 5e-14)
max_posterior_pdf = -Optim.minimum(result)
xpbdot_prior_max = Optim.minimizer(result)

# WE OBTAINED MU AND SIGMA FROM KDE PRIOR INFORMATION 
xpbdot_prior_kde_mu, _ = quadgk(x -> x * xpbdot_prior_kde(x), -5e-14, 5e-14)
xpbdot_prior_kde_sigma, _ = quadgk(x -> (x - xpbdot_prior_kde_mu)^2 * xpbdot_prior_kde(x), -5e-14, 5e-14) .^ (1/2)

#-------------------------------------------------------------------------------------
# NOW WE PREPARE TO BUILD POSTERIOR
# WE HAVE TRUE PRIOR AND NORMAL LIKELIHOOD -> WE CAN BUILD QUASI-NORMAL POSTERIOR
# WE APPROXIMATE LIKELIHOOD CALCULATING IT ON A RANGE OF XPBDOT
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# 0th APPROXIMATION BASED ON NORMAL PRIOR AND NORMAL LIKELIHOOD -> NORMAL POSTERIOR
#-------------------------------------------------------------------------------------

# XPBDOT RANGE CENTERED NEAR NORMAL POSTERIOR MAXIMUM
N_xpbdot_sigma = 5
xpbdot_posterior_sigma = (xpbdot_prior_sigma^-2 + xpbdot_sigma^-2)^(-1/2)
xpbdot_posterior_mu = xpbdot_posterior_sigma^2 * (xpbdot_prior_mu / xpbdot_prior_sigma^2 + xpbdot_mu / xpbdot_sigma^2)
xpbdot_posterior_range_inv = xpbdot_posterior_mu + N_xpbdot_sigma*xpbdot_posterior_sigma : -xpbdot_posterior_sigma : xpbdot_posterior_mu - (N_xpbdot_sigma + 0.1)*xpbdot_posterior_sigma
#--------------------------------------


#-------------------------------------------------------------------------------------
# 1st APPROXIMATION BASED ON TRUE PRIOR AND NORMAL LIKELIHOOD -> QUASI-NORMAL POSTERIOR
# A) XPBDOT RANGE CENTERED NEAR NORMAL LIKELIHOOD MAXIMUM
#-------------------------------------------------------------------------------------

# WE DO NOT NEED QUASI-NORMAL POSTERIOR

# BUILD RANGE
N_xpbdot_sigma = 5
xpbdot_range_inv = xpbdot_mu + N_xpbdot_sigma*xpbdot_sigma : -xpbdot_sigma : xpbdot_mu - (N_xpbdot_sigma + 0.1)*xpbdot_sigma
xpbdot_range = reverse(xpbdot_range_inv)
xpbdot_tempo_arr = collect(xpbdot_range)

global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 2*N_xpbdot_sigma+1,
    nits = 3 * ones(2*N_xpbdot_sigma+1),
    tparams_local = [[TP("XPBDOT", value, flag = 0)] for value in xpbdot_range_inv]
    )

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

chisqr_tempo_arr = reverse([results_global_iters.all_global_iterations[i].last_internal_iteration.result.chisqr for i in 1:global_iters_settings.iters])

#--------------------------------------
# RESULT for N_xpbdot_sigma = 5
chisqr_tempo_arr = [20949.67, 20936.43, 20926.39, 20919.47, 20915.51, 20914.19, 20914.87, 20916.53, 20918.38, 20920.15, 20921.86]
delta_chisqr_tempo_arr = chisqr_tempo_arr .- chisqr_min

#--------------------------------------
# version with sweep parameters
parameter_sweep_settings = ParameterSweepSettings(
    parameter_name = "XPBDOT",
    #values = collect(-3e-14:1e-15:3e-14)
    values = collect(xpbdot_range_inv)
)

results_sweep = run_tempo_parameter_sweep(basic_settings, parameter_sweep_settings)

chisqr_tempo_arr = [results_sweep[iter].last_internal_iteration.result.chisqr for iter in 1:length(parameter_sweep_settings.values)]
delta_chisqr_tempo_arr = chisqr_tempo_arr .- chisqr_min
#--------------------------------------


#-------------------------------------------------------------------------------------
# 1st APPROXIMATION BASED ON TRUE PRIOR AND NORMAL LIKELIHOOD -> QUASI-NORMAL POSTERIOR
# B) XPBDOT RANGE CENTERED NEAR APPROXIMATED POSTERIOR MAXIMUM
#-------------------------------------------------------------------------------------

#--------------------------------------


# BUILD APPROXIMATED QUASI-NORMAL POSTERIOR
product_prior_kde_normal_likelihood(xpbdot) = xpbdot_prior_kde(xpbdot) * xpbdot_normal_likelihood(xpbdot)
normalizing_constant_kde_likelihood, _ = quadgk(product_prior_kde_normal_likelihood, -5e-14, 5e-14)
quasi_normal_posterior(xpbdot) = product_prior_kde_normal_likelihood(xpbdot) / normalizing_constant_kde_likelihood

# WE FOUND QUASI-NORMAL POSTERIOR MAXIMUM
result = optimize(x -> -quasi_normal_posterior(x), -5e-14, 5e-14)
max_posterior_pdf = -Optim.minimum(result)
xpbdot_posterior_max = Optim.minimizer(result)

# WE OBTAINED MU AND SIGMA FOR QUASI-NORMAL POSTERIOR 
xpbdot_posterior_mu, _ = quadgk(x -> x * quasi_normal_posterior(x), -5e-14, 5e-14)
xpbdot_posterior_sigma, _ = quadgk(x -> (x - xpbdot_posterior_mu)^2 * quasi_normal_posterior(x), -5e-14, 5e-14) .^ (1/2)

# BUILD RANGE
N_xpbdot_posterior_sigma = 5
xpbdot_range_inv = xpbdot_posterior_mu + N_xpbdot_posterior_sigma*xpbdot_posterior_sigma : -xpbdot_posterior_sigma : xpbdot_posterior_mu - (N_xpbdot_posterior_sigma + 0.1)*xpbdot_posterior_sigma
xpbdot_range = reverse(xpbdot_range_inv)
xpbdot_tempo_arr = collect(xpbdot_range)

global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 2*N_xpbdot_posterior_sigma+1,
    nits = 3 * ones(2*N_xpbdot_posterior_sigma+1),
    tparams_local = [[TP("XPBDOT", value, flag = 0)] for value in xpbdot_range_inv]
    )

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

chisqr_tempo_arr = reverse([results_global_iters.all_global_iterations[i].last_internal_iteration.result.chisqr for i in 1:global_iters_settings.iters])


#--------------------------------------
# RESULT for N_xpbdot_sigma = 5 and xpbdot_prior .+= 0.0
chisqr_tempo_arr = [20917.04, 20915.82, 20914.96, 20914.42, 20914.2, 20914.24, 20914.5, 20914.94, 20915.49, 20916.11, 20916.77]
delta_chisqr_tempo_arr = chisqr_tempo_arr .- chisqr_min



# RESULT for N_xpbdot_sigma = 5 and xpbdot_prior .+= 1e-14
chisqr_tempo_arr = [20915.01, 20915.36, 20915.75, 20916.16, 20916.59, 20917.02, 20917.45, 20917.88, 20918.31, 20918.73, 20919.17]
delta_chisqr_tempo_arr = chisqr_tempo_arr .- chisqr_min




#--------------------------------------
# PLOT TRUE PRIOR, NORMAL LIKELIHOOD AND THEIR APPROXIMATED POSTERIOR
xpbdot_fine_arr = collect(-5e-14:5e-17:5e-14)
fig, ax = subplots()
plot(xpbdot_fine_arr, xpbdot_normal_likelihood.(xpbdot_fine_arr), label = "Normal likelihood")
#hist(xpbdot, bins=200, density=true)
plot(xpbdot_fine_arr, xpbdot_prior_kde(xpbdot_fine_arr), label = "True prior")
plot(xpbdot_fine_arr, quasi_normal_posterior.(xpbdot_fine_arr), label = "Quasi-normal posterior")
xlabel("XPBDOT, s/s")
ylabel("pdf")
legend()
#--------------------------------------


#--------------------------------------






#-------------------------------------------------------------------------------------
# PLOTS

N_points_arr = [3, 5, 7, 9, 11]
colors = ["b", "g", "r", "c", "m", "y", "k"]

fig, ax = subplots()

for (i, N_points) in enumerate(N_points_arr)
    N_points_side = div(N_points, 2)
    xpbdot_tempo_range_local = xpbdot_posterior_max - N_points_side * xpbdot_posterior_sigma : xpbdot_posterior_sigma : xpbdot_posterior_max + (N_points_side + 0.1) * xpbdot_posterior_sigma
    xpbdot_tempo_arr_local = collect(xpbdot_tempo_range_local)
    println(xpbdot_tempo_arr_local)
    N_full = length(delta_chisqr_tempo_arr)
    i_left, i_right = 1 + div(N_full,2) - N_points_side, 1 + div(N_full,2) + N_points_side
    delta_chisqr_spline = cubic_spline_interpolation(xpbdot_tempo_range_local, delta_chisqr_tempo_arr[i_left:i_right], extrapolation_bc = Line())

    likelihood_spline(xpbdot) = exp(-0.5 * delta_chisqr_spline(xpbdot))
    normalizing_constant_likelihood, _ = quadgk(likelihood_spline, -5e-14, 5e-14)
    likelihood(xpbdot) = likelihood_spline(xpbdot) / normalizing_constant_likelihood
    plot(xpbdot_tempo_arr_local, likelihood.(xpbdot_tempo_arr_local), "d", color=colors[i])
    plot(xpbdot_fine_arr, likelihood.(xpbdot_fine_arr), label = "Splines likelihood $(N_points)", color=colors[i])

end

plot(xpbdot_fine_arr, xpbdot_normal_likelihood.(xpbdot_fine_arr), label = "Normal likelihood", color = "k")

xlabel("XPBDOT, s/s")
ylabel("pdf")
legend()

#--------------------------------------

N_points_arr = [3, 5, 7, 9, 11]
colors = ["b", "g", "r", "c", "m", "y", "k"]

fig, ax = subplots()

for (i, N_points) in enumerate(N_points_arr)
    N_points_side = div(N_points, 2)
    xpbdot_tempo_range_local = xpbdot_posterior_max - N_points_side * xpbdot_posterior_sigma : xpbdot_posterior_sigma : xpbdot_posterior_max + (N_points_side + 0.1) * xpbdot_posterior_sigma
    xpbdot_tempo_arr_local = collect(xpbdot_tempo_range_local)
    println(xpbdot_tempo_arr_local)
    N_full = length(delta_chisqr_tempo_arr)
    i_left, i_right = 1 + div(N_full,2) - N_points_side, 1 + div(N_full,2) + N_points_side
    delta_chisqr_spline = cubic_spline_interpolation(xpbdot_tempo_range_local, delta_chisqr_tempo_arr[i_left:i_right], extrapolation_bc = Line())

    result = optimize(x -> delta_chisqr_spline(x), -5e-14, 5e-14)
    delta_chisqr_spline_min = Optim.minimum(result)

    plot(xpbdot_tempo_arr_local, delta_chisqr_spline.(xpbdot_tempo_arr_local), "d", color=colors[i])
    plot(xpbdot_fine_arr, delta_chisqr_spline.(xpbdot_fine_arr), label = "Splines likelihood $(N_points)", color=colors[i])

end

plot(xpbdot_fine_arr, delta_chisqr_normal_likelihood.(xpbdot_fine_arr), label = "Normal likelihood", color = "k")

xlabel("XPBDOT, s/s")
ylabel(L"\Delta \chi^2")
legend()



#--------------------------------------

N_points_arr = [3, 5, 7, 9, 11]
colors = ["b", "g", "r", "c", "m", "y", "k"]

fig, ax = subplots()

for (i, N_points) in enumerate(N_points_arr)
    N_points_side = div(N_points, 2)
    xpbdot_tempo_range_local = xpbdot_posterior_max - N_points_side * xpbdot_posterior_sigma : xpbdot_posterior_sigma : xpbdot_posterior_max + (N_points_side + 0.1) * xpbdot_posterior_sigma
    xpbdot_tempo_arr_local = collect(xpbdot_tempo_range_local)
    println(xpbdot_tempo_arr_local)
    N_full = length(delta_chisqr_tempo_arr)
    i_left, i_right = 1 + div(N_full,2) - N_points_side, 1 + div(N_full,2) + N_points_side
    delta_chisqr_spline = cubic_spline_interpolation(xpbdot_tempo_range_local, delta_chisqr_tempo_arr[i_left:i_right], extrapolation_bc = Line())

    result = optimize(x -> delta_chisqr_spline(x), -5e-14, 5e-14)
    delta_chisqr_spline_min = Optim.minimum(result)

    temp_fun = xpbdot -> exp(-0.5 * delta_chisqr_spline(xpbdot))
    normalizing_constant_likelihood, _ = quadgk(temp_fun, -5e-14, 5e-14)
    likelihood_spline(xpbdot) = temp_fun(xpbdot) / normalizing_constant_likelihood

    product_prior_kde_likelihood_spline(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_spline(xpbdot)
    normalizing_kde_likelihood_spline, _ = quadgk(product_prior_kde_likelihood_spline, -5e-14, 5e-14)
    posterior_spline(xpbdot) = product_prior_kde_likelihood_spline(xpbdot) / normalizing_kde_likelihood_spline

    plot(xpbdot_tempo_arr_local, posterior_spline.(xpbdot_tempo_arr_local), "d", color=colors[i])
    plot(xpbdot_fine_arr, posterior_spline.(xpbdot_fine_arr), label = "Splines posterior $(N_points)", color=colors[i])

end

plot(xpbdot_fine_arr, quasi_normal_posterior.(xpbdot_fine_arr), label = "Quasi-normal posterior", color = "k")

xlabel("XPBDOT, s/s")
ylabel("pdf")
legend()



#--------------------------------------
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

function perform_marginalization(prior::Function, likelihood::Function)
    normalizing_constant_likelihood, _ = quadgk(likelihood, -5e-14, 5e-14)
    P, _ = quadgk(x -> prior(x) * likelihood(x) / normalizing_constant_likelihood, -5e-14, 5e-14)

    result = optimize(x -> - likelihood(x), -5e-14, 5e-14)
    xpbdot_likelihood_max = Optim.minimizer(result)

    P_likelihood_max = likelihood(xpbdot_likelihood_max) / normalizing_constant_likelihood
    delta_chisqr_marginalized = -2 * log(P) + 2 * log(P_likelihood_max)
    return delta_chisqr_marginalized
end


N_points_arr = [3, 5, 7, 9, 11]
colors = ["b", "g", "r", "c", "m", "y", "k"]

fig, ax = subplots()

for (i, N_points) in enumerate(N_points_arr)
    N_points_side = div(N_points, 2)
    xpbdot_tempo_range_local = xpbdot_posterior_max - N_points_side * xpbdot_posterior_sigma : xpbdot_posterior_sigma : xpbdot_posterior_max + (N_points_side + 0.1) * xpbdot_posterior_sigma
    xpbdot_tempo_arr_local = collect(xpbdot_tempo_range_local)
    # println(xpbdot_tempo_arr_local)
    N_full = length(delta_chisqr_tempo_arr)
    i_left, i_right = 1 + div(N_full,2) - N_points_side, 1 + div(N_full,2) + N_points_side
    delta_chisqr_spline = cubic_spline_interpolation(xpbdot_tempo_range_local, delta_chisqr_tempo_arr[i_left:i_right], extrapolation_bc = Line())

    result = optimize(x -> delta_chisqr_spline(x), -5e-14, 5e-14)
    delta_chisqr_spline_min = Optim.minimum(result)
    xpbdot_delta_chisqr_spline_min = Optim.minimizer(result)

    temp_fun = xpbdot -> exp(-0.5 * delta_chisqr_spline(xpbdot))
    normalizing_constant_likelihood, _ = quadgk(temp_fun, -5e-14, 5e-14)
    likelihood_spline(xpbdot) = temp_fun(xpbdot) / normalizing_constant_likelihood

    product_prior_kde_likelihood_spline(xpbdot) = xpbdot_prior_kde(xpbdot) * likelihood_spline(xpbdot)
    normalizing_kde_likelihood_spline, _ = quadgk(product_prior_kde_likelihood_spline, -5e-14, 5e-14)
    posterior_spline(xpbdot) = product_prior_kde_likelihood_spline(xpbdot) / normalizing_kde_likelihood_spline

    delta_chisqr_posterior_fine_arr, max_posterior_xpbdot, posterior_correction = get_delta_chisqr_posterior(xpbdot_fine_arr, xpbdot_prior_kde, posterior_spline)

    full_correction = delta_chisqr_spline(max_posterior_xpbdot) + posterior_correction

    plot(max_posterior_xpbdot, full_correction, "d", color=colors[i])
    plot(xpbdot_fine_arr, delta_chisqr_spline.(xpbdot_fine_arr) + delta_chisqr_posterior_fine_arr, label = "Splines $(N_points) " * string(full_correction), color=colors[i])

    delta_chisqr_marginalized = perform_marginalization(xpbdot_prior_kde, posterior_spline)
    println("Splines $(N_points). Correction: max posterior = $(full_correction), marginalized = $(delta_chisqr_marginalized)")
end

delta_chisqr_posterior_fine_arr, max_posterior_xpbdot, posterior_correction = get_delta_chisqr_posterior(xpbdot_fine_arr, xpbdot_prior_kde, quasi_normal_posterior)
full_correction = delta_chisqr_normal_likelihood.(max_posterior_xpbdot) + posterior_correction

plot(xpbdot_fine_arr, delta_chisqr_normal_likelihood.(xpbdot_fine_arr) + delta_chisqr_posterior_fine_arr, label = "Quasi-normal posterior " * string(full_correction), color="k")


xlabel("XPBDOT, s/s")
ylabel(L"\Delta \chi^2 \mathrm{, \, posterior}")
legend()

#--------------------------------------
