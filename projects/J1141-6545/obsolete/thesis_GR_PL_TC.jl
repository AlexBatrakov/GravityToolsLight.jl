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

#-------------------------------------------------------------------------------------

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551

N_window = 5
N_iters = N_window * 46 + 1
time_arr = collect(LinRange(time_start, time_finish, N_iters))

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_PL_TC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PL_TC.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals -epoch center",
    tparams = [TP("NITS", 20), TP("GAIN_ONE_AFTER_NITS", 10)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = N_iters - N_window,
    nits = [5],
    gain = [1],
    tparams_local = [[TP("START", time_arr[i], flag=1), TP("FINISH", time_arr[i + N_window], flag=1)] for i in 1:N_iters-N_window],
    # flags = [""],
    # fit_EFACs_EQUADs = [false],
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
                if param_name == :offset
                    push!(param_values, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.offset[1])
                    push!(param_errors, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.offset[2])
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

params_values = Dict{Symbol, Vector{Float64}}()
params_errors = Dict{Symbol, Vector{Float64}}()

# param_names = [:PEPOCH, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :M2, :MTOT, :I, :offset, :AR, :OMDOT, :PBDOT, :GAMMA]
# param_names = [:PEPOCH, :F0, :F1, :F2, :offset]
param_names = [:PEPOCH, :F0, :F1, :offset]

for param_name in param_names
    params_values[param_name], params_errors[param_name] = extract_parameter(results_global_iters, param_name)
end

for param_name in [:offset, :F0, :F1]
    epoch = 55825.0
    t = params_values[:PEPOCH]
    values = params_values[param_name]
    errors = params_errors[param_name]
    weights = errors .^ (-2.0)
    mean_weighted_param_value = sum(weights .* values) / sum(weights)
    sigma_mean = sqrt(1.0 / sum(weights))

    p = Polynomials.fit(t, values, 1, weights = weights)  # Полином степени 1
    trend = p.(t)

    fig, ax = subplots()
    errorbar(t, values, yerr = errors, fmt=".")
    plot(t, trend, label="linear trend")
    # Форматируем строку с учётом ошибок
    params_text = string(
        @sprintf("Weighted mean: %.8e ± %.8e", mean_weighted_param_value, sigma_mean)
    )
    # Координаты для текста на графике
    x_pos = 0.05  # Процентное смещение по оси X (0.05 = 5% от левого края)
    y_pos = 0.05  # Процентное смещение по оси Y (0.95 = 95% от нижнего края)
    # Добавляем текст с плашкой
    text(
        x_pos, y_pos, params_text, transform=gca().transAxes,
        fontsize=10,
        bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    )
    xlabel("MJD")
    ylabel("$(param_name)")
    grid()
    legend()
    tight_layout()
    show()
end

#---------------------------------------------------------------------------------------------------------------

param_name = :F0
T_EP = 55825.0
t = params_values[:PEPOCH] .- T_EP
y = params_values[param_name]
y_err = params_errors[param_name]
y_err_median = median(y_err)
ind = y_err .< Inf * y_err_median

t = t[ind]
y = y[ind]
y_err = y_err[ind]
y_weights = (y_err) .^ (-2.0)
# weights = 1.0 ./ abs.(y_err)
poly_deg = 2

pol_F = Polynomials.fit(t, y, poly_deg, weights = y_weights)  # Полином степени 2

# Предсказанный тренд
trend = pol_F.(t)

# Данные без тренда
detrended_y = y .- trend
detrended_y_abs_median = median(abs.(detrended_y))

fig, ax = subplots()
errorbar(t .+ T_EP, detrended_y, yerr = y_err, fmt=".")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

#---------------------------------------------------------------------------------------------------------------


using LombScargle
pgram = lombscargle(t, detrended_y, y_err)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")
xlim(left=(time_finish - time_start) / (N_iters - 1))


#---------------------------------------------------------------------------------------------------------------

function v_x(t, T0, Pb, e, x, om)
    Pb, e, x = abs.((Pb, e, x))
    e = min(e, 1.0)
    M = 2.0 * pi / Pb * (t - T0)
    U_func = U -> U - e * sin(U) - M
    U = find_zero(U_func, M)

    sin_AT = sqrt(1.0 - e^2) * sin(U) / (1.0 - e * cos(U))
    cos_AT = (cos(U) - e) / (1.0 - e * cos(U))
    AT = atan(sin_AT, cos_AT)

    v_x = 2.0 * pi / (Pb * 86400.0) * x / sqrt(1 - e^2) * (cos(om + AT) + e * cos(om))


    # dU_dt = 2.0 * pi / Pb / (1 - e * cos(U)) 
    # dRoe = x * sin(om) * (cos(U) - e) + x * sqrt(1 - e^2) * cos(om) * sin(U)
    # v_x = x * sin(om) * sin(U) * dU_dt - x * sqrt(1 - e^2) * cos(om) * cos(U) * dU_dt
    return v_x
end

planet_F_model(t, p) = (1.0 .- v_x.(t, p[1] - T_EP, p[2], p[3], p[4], p[5] * pi/180.0)) .* (p[6] .+ p[7] .* t .+ p[8] .* t .^ 2)

T0_est = 55650.0
Pb_est = 1284.5
e_est  = 0.2
x_est  = 0.0001
om_est = 0

p0 = [T0_est, Pb_est, e_est, x_est, om_est, pol_F[0], pol_F[1], pol_F[2]]

est_trend = pol_F.(t)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
est_trend_fine = pol_F.(t_fine)
est_signal_fine = planet_F_model(t_fine, p0)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- est_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, est_signal_fine .- est_trend_fine, label="Planet est")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))


fit_result = curve_fit(planet_F_model, t, y, y_weights, p0)

T0, Pb, e, x, omega0, F0, F1, F2 = fit_result.param
dT0, dPb, de, dx, domega0, dF0, dF1, dF2 = standard_errors(fit_result)

predicted_trend = F0 .+ F1 .* t .+ F2 .* t .^2
predicted_signal = planet_F_model(t, fit_result.param)

predicted_trend_fine = F0 .+ F1 .* t_fine .+ F2 .* t_fine .^2
predicted_signal_fine = planet_F_model(t_fine, fit_result.param)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fittes planet")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))
seconds_per_day = 86400.0

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("PL_T0  = %.4f ± %.4f", T0, dT0), "\n",
    @sprintf("PL_PB  = %.2f ± %.2f", Pb, dPb), "\n",
    @sprintf("PL_ECC = %.5f ± %.5f",  e, de), "\n",
    @sprintf("PL_X   = %.5f ± %.5f",  x, dx), "\n",
    @sprintf("PL_OM  = %.2f ± %.2f",  omega0, domega0), "\n",
    @sprintf("F0     = %.11f ± %.2e", F0, dF0), "\n",
    @sprintf("F1     = %.8e ± %.2e", F1 / seconds_per_day, dF1 / seconds_per_day), "\n",
    @sprintf("F2     = %.5e ± %.2e", F2 / seconds_per_day^2, dF2 / seconds_per_day^2)
)

# Координаты для текста на графике
x_pos = 0.05  # Процентное смещение по оси X (0.05 = 5% от левого края)
y_pos = 0.05  # Процентное смещение по оси Y (0.95 = 95% от нижнего края)
# Добавляем текст с плашкой
text(
    x_pos, y_pos, params_text, transform=gca().transAxes,
    fontsize=10,
    bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
)

MJD_glitch = 54272.697451999999998
axvline(x=MJD_glitch, color="red", linestyle="--", linewidth=1, label="Glitch epoch")

# Оформление
xlabel("MJD")
ylabel(L"\Delta F\ (s^{-1})")
legend(loc="upper right")
title(L"\Delta F = F - (F_0 + F_1 (T - T_0) + \frac{1}{2} F_2 (T - T_0)^2)")
grid()
tight_layout()
show()

#---------------------------------------------------------------------------------------------------------------

residual_signal = y .- predicted_signal
fig, ax = subplots()
errorbar(t, residual_signal, yerr = y_err, fmt=".")

#---------------------------------------------------------------------------------------------------------------

function v_x(t, T0, Pb, e, x, om)
    Pb, e, x = abs.((Pb, e, x))
    e = min(e, 1.0)
    M = 2.0 * pi / Pb * (t - T0)
    U_func = U -> U - e * sin(U) - M
    U = find_zero(U_func, M)

    sin_AT = sqrt(1.0 - e^2) * sin(U) / (1.0 - e * cos(U))
    cos_AT = (cos(U) - e) / (1.0 - e * cos(U))
    AT = atan(sin_AT, cos_AT)

    v_x = 2.0 * pi / (Pb * 86400.0) * x / sqrt(1 - e^2) * (cos(om + AT) + e * cos(om))


    # dU_dt = 2.0 * pi / Pb / (1 - e * cos(U)) 
    # dRoe = x * sin(om) * (cos(U) - e) + x * sqrt(1 - e^2) * cos(om) * sin(U)
    # v_x = x * sin(om) * sin(U) * dU_dt - x * sqrt(1 - e^2) * cos(om) * cos(U) * dU_dt
    return v_x
end

planet_F_model(t, p) = (1.0 .- v_x.(t, p[1] - T_EP, p[2], p[3], p[4], p[5] * pi/180.0) .- v_x.(t, p[6] - T_EP, p[7], p[8], p[9], p[10] * pi/180.0)) .* (p[11] .+ p[12] .* t .+ p[13] .* t .^ 2)

T0_est1 = 55260.0
Pb_est1 = 750.0
e_est1  = 0.
x_est1  = 0.0004
om_est1 = 0

T0_est2 = 55260.0
Pb_est2 = 1284.5
e_est2  = 0.
x_est2  = 0.0004
om_est2 = 0

p0 = [T0_est1, Pb_est1, e_est1, x_est1, om_est1, T0_est2, Pb_est2, e_est2, x_est2, om_est2, pol_F[0], pol_F[1], pol_F[2]]

est_trend = pol_F.(t)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
est_trend_fine = pol_F.(t_fine)
est_signal_fine = planet_F_model(t_fine, p0)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- est_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, est_signal_fine .- est_trend_fine, label="Planet est")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))


fit_result = curve_fit(planet_F_model, t, y, y_weights, p0)

T01, Pb1, e1, x1, omega01, T02, Pb2, e2, x2, omega02, F0, F1, F2 = fit_result.param
dT0, dPb, de, dx, domega0, dF0, dF1, dF2 = standard_errors(fit_result)

predicted_trend = F0 .+ F1 .* t .+ F2 .* t .^2
predicted_signal = planet_F_model(t, fit_result.param)

predicted_trend_fine = F0 .+ F1 .* t_fine .+ F2 .* t_fine .^2
predicted_signal_fine = planet_F_model(t_fine, fit_result.param)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fittes planet")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))
seconds_per_day = 86400.0

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("PL_T0  = %.4f ± %.4f", T0, dT0), "\n",
    @sprintf("PL_PB  = %.2f ± %.2f", Pb, dPb), "\n",
    @sprintf("PL_ECC = %.5f ± %.5f",  e, de), "\n",
    @sprintf("PL_X   = %.5f ± %.5f",  x, dx), "\n",
    @sprintf("PL_OM  = %.2f ± %.2f",  omega0, domega0), "\n",
    @sprintf("F0     = %.11f ± %.2e", F0, dF0), "\n",
    @sprintf("F1     = %.8e ± %.2e", F1 / seconds_per_day, dF1 / seconds_per_day), "\n",
    @sprintf("F2     = %.5e ± %.2e", F2 / seconds_per_day^2, dF2 / seconds_per_day^2)
)

# Координаты для текста на графике
x_pos = 0.05  # Процентное смещение по оси X (0.05 = 5% от левого края)
y_pos = 0.05  # Процентное смещение по оси Y (0.95 = 95% от нижнего края)
# Добавляем текст с плашкой
text(
    x_pos, y_pos, params_text, transform=gca().transAxes,
    fontsize=10,
    bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
)

MJD_glitch = 54272.697451999999998
axvline(x=MJD_glitch, color="red", linestyle="--", linewidth=1, label="Glitch epoch")

# Оформление
xlabel("MJD")
ylabel(L"\Delta F\ (s^{-1})")
legend(loc="upper right")
title(L"\Delta F = F - (F_0 + F_1 (T - T_0) + \frac{1}{2} F_2 (T - T_0)^2)")
grid()
tight_layout()
show()