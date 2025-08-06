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
pygui(true)

#-------------------------------------------------------------------------------------

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551

N_iters = 47
time_arr = collect(LinRange(time_start, time_finish, N_iters))

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_TC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_TC.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals -epoch center",
    tparams = [TP("NITS", 6)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = N_iters - 1,
    nits = [6],
    gain = [1],
    tparams_local = [[TP("START", time_arr[i], flag=1), TP("FINISH", time_arr[i+1], flag=1)] for i in 1:N_iters-1],
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
param_names = [:PEPOCH, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :I, :offset]

for param_name in param_names
    params_values[param_name], params_errors[param_name] = extract_parameter(results_global_iters, param_name)
end

params_values[:I] = [I < 90.0 ? I : 180.0 - I for I in params_values[:I]]

# params_values, params_errors = extract_parameters(results_global_iters, [:PEPOCH, :F0, :F1, :PB, :A1, :OM, :ECC, :I])

# PEPOCH_values, _ = extract_parameter(results_global_iters, :PEPOCH)
# F0_values, F0_errors = extract_parameter(results_global_iters, :F0)
# F1_values, F1_errors = extract_parameter(results_global_iters, :F1)
# PB_values, PB_errors = extract_parameter(results_global_iters, :PB)
# A1_values, A1_errors = extract_parameter(results_global_iters, :A1)
# OM_values, OM_errors = extract_parameter(results_global_iters, :OM)
# ECC_values, ECC_errors = extract_parameter(results_global_iters, :ECC)

for param_name in [:F0, :F1, :F2]
    epoch = 55825.0
    t = params_values[:PEPOCH]
    values = params_values[param_name]
    errors = params_errors[param_name]
    weights = errors .^ (-2.0)
    mean_weighted_param_value = sum(weights .* values) / sum(weights)
    sigma_mean = sqrt(1.0 / sum(weights))

    p = Polynomials.fit(t, values, 1, weights = weights)  # Полином степени 1
    trend = p.(t)

    # linear_model(t, p) = p[1] .* (t .- epoch) .+ p[2]
    # # Начальные параметры: [A, omega, T0, offset]
    # p0 = [0, mean_weighted_param_value]

    # # Подгонка
    # fit_result = curve_fit(linear_model, t, values, weights, p0)

    # # Коэффициенты
    # inclination, offset = fit_result.param
    # dinclination, doffset = standard_errors(fit_result)

    # trend = linear_model(t, fit_result.param)

    fig, ax = subplots()
    errorbar(t, values, yerr = errors, fmt=".")
    plot(t, trend, label="linear trend")
    # Форматируем строку с учётом ошибок
    params_text = string(
        @sprintf("Weighted mean: %.8e ± %.8e", mean_weighted_param_value, sigma_mean)
        # @sprintf("Linear mean: %.8e ± %.8e", offset, doffset), "\n",
        # @sprintf("Linear inclination: %.8e ± %.8e", inclination, dinclination)
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





# using DataFrames
# t = params_values[:PEPOCH] .- time_start
# y = params_values[:F0]
# data = DataFrame(t=t, y=y)

# # Линейная модель
# model = lm(@formula(y ~ t + t^2), data)

# # Коэффициенты
# coeffs = coef(model)

# # Предсказанный тренд
# trend = predict(model)

# # Данные без тренда
# detrended_y = y .- trend

# model = lm(@formula(y ~ t + t*t), data)

#-----------------------------------------------
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
weights = (y_err) .^ (-2.0)
poly_deg = 2

p = Polynomials.fit(t, y, poly_deg, weights = weights)  # Полином степени 2

# Предсказанный тренд
trend = p.(t)

# Данные без тренда
detrended_y = y .- trend
detrended_y_abs_median = median(abs.(detrended_y))

fig, ax = subplots()
errorbar(t .+ T_EP, detrended_y, yerr = y_err, fmt=".")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

#-----------------------------------------------


using LombScargle
pgram = lombscargle(t, detrended_y, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")
xlim(left=(time_finish - time_start) / (N_iters - 1))
# axvline(x=(time_finish - time_start) / (N_iters - 1), color="gray", linestyle="--")
P1_estimate = 3400.0
P2_estimate = 860.0

#-----------------------------------------------

using LsqFit
using Printf
# Нелинейная модель
one_harminic_model(t, p) = - p[1] .* cos.(p[2] .* (t .- p[3])) .+ p[4] .+ p[5] .* t .+ p[6] .* t .^ 2

# Начальные параметры: [A, omega, T0, offset]
p0 = [detrended_y_abs_median, 2*pi/P1_estimate, 55500.0 - T_EP, p[0], p[1], p[2]]

# Подгонка
fit_result = curve_fit(one_harminic_model, t, y, weights, p0)

# Коэффициенты
A, omega, T0, a, b, c = fit_result.param
dA, domega, dT0, da, db, dc = standard_errors(fit_result)
P = 2*pi/omega
dP = abs(- 2*pi/omega^2 * domega)
println("Амплитуда косинуса (A): $A ($dA)")
println("Угловая частота (omega): $omega ($domega)")
println("Период (P): $(P) ($dP)")
println("Время нулевой фазы (T0): $T0 ($dT0)")
println("a: $a ($da)")
println("b: $b ($db)")
println("c: $c ($dc)")

# Предсказанные значения
predicted_trend = a .+ b .* t .+ c .* t .^2
predicted_signal = one_harminic_model(t, fit_result.param)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
predicted_trend_fine = a .+ b .* t_fine .+ c .* t_fine .^2
predicted_signal_fine = one_harminic_model(t_fine, fit_result.param)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fitted harmonic")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

seconds_per_day = 86400.0

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("A = %.4e ± %.2e", A, dA), "\n",
    @sprintf("P = %.2f ± %.2f", P, dP), "\n",
    @sprintf("T0 = %.2f ± %.2f", T0 + T_EP, dT0), "\n",
    @sprintf("F0 = %.11f ± %.2e", a, da), "\n",
    @sprintf("F1 = %.8e ± %.2e", b / seconds_per_day, db / seconds_per_day), "\n",
    @sprintf("F2 = %.5e ± %.2e", c / seconds_per_day^2, dc / seconds_per_day^2)
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

#-----------------------------------------------

# pgram = lombscargle(t, y .- predicted_trend, y_err, normalization=:psd)
# fig, ax = subplots()
# plot(periodpower(pgram)...)
# yscale("log")
# xscale("log")

#-----------------------------------------------

residual_signal = y .- predicted_signal
fig, ax = subplots()
errorbar(t, residual_signal, yerr = y_err, fmt=".")

using LombScargle
pgram = lombscargle(t, residual_signal, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")


#-----------------------------------------------
# # Нелинейная модель
# two_harminic_model(t, p) = p[1] .* cos.(p[2] .* (t .- p[3])) .+ p[4] .* cos.(p[5] .* (t .- p[6])) .+ p[7] .+ p[8] .* t .+ p[9] .* t .^ 2

# # Начальные параметры: [A, omega, T0, offset]
# p0 = [A, 2*pi/P, T0, 0.1 * A, 2*pi/P2_estimate, 55500.0 - T_EP, p[0], p[1], p[2]]

# # Подгонка
# fit_result = curve_fit(two_harminic_model, t, y, weights, p0)

# # Коэффициенты
# A1, omega1, T01, A2, omega2, T02, a, b, c = fit_result.param
# dA1, domega1, dT01, dA2, domega2, dT02, da, db, dc = standard_errors(fit_result)
# P1 = 2*pi/omega1
# dP1 = abs(- 2*pi/omega1^2 * domega1)
# P2 = 2*pi/omega2
# dP2 = abs(- 2*pi/omega2^2 * domega2)
# println("1: Амплитуда косинуса (A): $A1 ($dA1)")
# println("1: Угловая частота (omega): $omega1 ($domega1)")
# println("1: Период (P): $(P1) ($dP1)")
# println("1: Время нулевой фазы (T0): $T02 ($dT02)")
# println("2: Амплитуда косинуса (A): $A2 ($dA2)")
# println("2: Угловая частота (omega): $omega2 ($domega2)")
# println("2: Период (P): $(P2) ($dP2)")
# println("2: Время нулевой фазы (T0): $T02 ($dT02)")
# println("a: $a ($da)")
# println("b: $b ($db)")
# println("c: $c ($dc)")

# # Предсказанные значения
# predicted_trend = a .+ b .* t .+ c .* t .^2
# predicted_signal = two_harminic_model(t, fit_result.param)

# t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
# predicted_trend_fine = a .+ b .* t_fine .+ c .* t_fine .^2
# predicted_signal_fine = two_harminic_model(t_fine, fit_result.param)

# fig, ax = subplots()
# errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
# plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fitted harmonic")
# ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

# seconds_per_day = 86400.0

# # Форматируем строку с учётом ошибок
# params_text = string(
#     @sprintf("A1 = %.4e ± %.2e", A1, dA1), "\n",
#     @sprintf("P1 = %.2f ± %.2f", P1, dP1), "\n",
#     @sprintf("T01 = %.2f ± %.2f", T01 + T_EP, dT01), "\n",
#     @sprintf("A2 = %.4e ± %.2e", A2, dA2), "\n",
#     @sprintf("P2 = %.2f ± %.2f", P2, dP2), "\n",
#     @sprintf("T02 = %.2f ± %.2f", T02 + T_EP, dT02), "\n",
#     @sprintf("F0 = %.11f ± %.2e", a, da), "\n",
#     @sprintf("F1 = %.8e ± %.2e", b / seconds_per_day, db / seconds_per_day), "\n",
#     @sprintf("F2 = %.5e ± %.2e", c / seconds_per_day^2, dc / seconds_per_day^2)
# )

# # Координаты для текста на графике
# x_pos = 0.05  # Процентное смещение по оси X (0.05 = 5% от левого края)
# y_pos = 0.05  # Процентное смещение по оси Y (0.95 = 95% от нижнего края)
# # Добавляем текст с плашкой
# text(
#     x_pos, y_pos, params_text, transform=gca().transAxes,
#     fontsize=10,
#     bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
# )

# MJD_glitch = 54272.697451999999998
# axvline(x=MJD_glitch, color="red", linestyle="--", linewidth=1, label="Glitch epoch")

# # Оформление
# xlabel("MJD")
# ylabel(L"\Delta F\ (s^{-1})")
# legend(loc = "upper right")
# title(L"\Delta F = F - (F_0 + F_1 (T - T_0) + \frac{1}{2} F_2 (T - T_0)^2)")
# grid()
# tight_layout()
# show()

# #-----------------------------------------------

# residual_signal = y .- predicted_signal
# fig, ax = subplots()
# errorbar(t, residual_signal, yerr = y_err, fmt=".")

# using LombScargle
# pgram = lombscargle(t, residual_signal, y_err, normalization=:psd)
# fig, ax = subplots()
# plot(periodpower(pgram)...)
# yscale("log")
# xscale("log")

#-----------------------------------------------
param_name = :F1
t = params_values[:PEPOCH] .- T_EP
y = params_values[param_name]
y_err = params_errors[param_name]
y_err_median = median(y_err)
ind = y_err .< Inf * y_err_median

t = t[ind]
y = y[ind]
y_err = y_err[ind]
weights = (y_err) .^ (-2.0)
poly_deg = 1

p = Polynomials.fit(t, y, poly_deg, weights = weights)  # Полином степени 2

# Предсказанный тренд
trend = p.(t)

# Данные без тренда
detrended_y = y .- trend
detrended_y_abs_median = median(abs.(detrended_y))

fig, ax = subplots()
errorbar(t, detrended_y, yerr = y_err, fmt=".")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

#-----------------------------------------------


using LombScargle
pgram = lombscargle(t, detrended_y, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")

P1_estimate = 3450.0
P2_estimate = 1970.0

#-----------------------------------------------

using LsqFit
using Printf
# Нелинейная модель
one_harminic_model(t, p) = p[1] .* cos.(p[2] .* (t .- p[3])) .+ p[4] .+ p[5] .* t

# Начальные параметры: [A, omega, T0, offset]
p0 = [detrended_y_abs_median, 2*pi/P1_estimate, 55500.0, p[0], p[1]]

# Подгонка
fit_result = curve_fit(one_harminic_model, t, y, weights, p0)

# Коэффициенты
A, omega, T0, a, b = fit_result.param
dA, domega, dT0, da, db = standard_errors(fit_result)
P = 2*pi/omega
dP = abs(- 2*pi/omega^2 * domega)
println("Амплитуда косинуса (A): $A ($dA)")
println("Угловая частота (omega): $omega ($domega)")
println("Период (P): $(P) ($dP)")
println("Время нулевой фазы (T0): $T0 ($dT0)")
println("a: $a ($da)")
println("b: $b ($db)")

# Предсказанные значения
predicted_trend = a .+ b .* t
predicted_signal = one_harminic_model(t, fit_result.param)

t_fine = collect(LinRange(time_start, time_finish, 1000))
predicted_trend_fine = a .+ b .* t_fine
predicted_signal_fine = one_harminic_model(t_fine, fit_result.param)

fig, ax = subplots()
errorbar(t, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine, predicted_signal_fine .- predicted_trend_fine, label="Подогнанная синусоида")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

seconds_per_day = 86400.0

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("A = %.2e ± %.2e", A, dA), "\n",
    @sprintf("P = %.2f ± %.2f", P, dP), "\n",
    @sprintf("T0 = %.2f ± %.2f", T0, dT0), "\n",
    @sprintf("F1 = %.2e ± %.2e", a, da), "\n",
    @sprintf("F2 = %.2e ± %.2e", b / seconds_per_day, db / seconds_per_day)
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

# Оформление
xlabel("Time (t)")
ylabel("Signal (y)")
legend()
title("Подгонка синусоиды с LsqFit")
grid()
tight_layout()
show()

#-----------------------------------------------

pgram = lombscargle(t, y .- predicted_trend, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")

#-----------------------------------------------

residual_signal = y .- predicted_signal
fig, ax = subplots()
errorbar(t, residual_signal, yerr = y_err, fmt=".")

using LombScargle
pgram = lombscargle(t, residual_signal, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
weights = (y_err) .^ (-2.0)
# weights = 1.0 ./ abs.(y_err)
poly_deg = 2

p = Polynomials.fit(t, y, poly_deg, weights = weights)  # Полином степени 2

# Предсказанный тренд
trend = p.(t)

# Данные без тренда
detrended_y = y .- trend
detrended_y_abs_median = median(abs.(detrended_y))

fig, ax = subplots()
errorbar(t .+ T_EP, detrended_y, yerr = y_err, fmt=".")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

#-----------------------------------------------


using LombScargle
pgram = lombscargle(t, detrended_y, y_err, normalization=:psd)
fig, ax = subplots()
plot(periodpower(pgram)...)
yscale("log")
xscale("log")
#-----------------------------------------------


function v_x(T, T0, Pb, e, x, omega0)
    Pb, e, x = abs.((Pb, e, x))
    e = min(e, 1.0)
    M = 2 * pi / Pb * (T - T0)
    U_func = U -> U - e * sin(U) - M
    U = find_zero(U_func, M)
    omega = omega0 * pi / 180.0
    dU_dt = 2 * pi / (Pb * 86400.0) / (1 - e * cos(U)) 
    dRoe = x * sin(omega) * (cos(U) - e) + x * sqrt(1 - e^2) * cos(omega) * sin(U)
    v_x = x * sin(omega) * sin(U) * dU_dt - x * sqrt(1 - e^2) * cos(omega) * cos(U) * dU_dt
    return v_x
end

using LsqFit
using Printf
# Нелинейная модель
planet_model(t, p) = (1.0 .+ v_x.(t, p[1], p[2], p[3], p[4], p[5])) .* (p[6] .+ p[7] .* t .+ p[8] .* t .^ 2)

# Начальные параметры: [A, omega, T0, offset]
p0 = [55000 - T_EP, 3400.0, 0.2, 0.008, 0.0, p[0], p[1], p[2]]

# Подгонка
fit_result = curve_fit(planet_model, t, y, weights, p0)

# function planet_model_objective(p, t=t, y=y, err=y_err)
#     T0, Pb, e, x, omega0, F0, F1, F2 = p
#     res = y .- planet_model(t, p)
#     return sqrt(sum((res ./ y_err) .^ 2 ./ length(t)))
# end

# optim_result = optimize(planet_model_objective, p0)
# optim_result = optimize(planet_model_objective, optim_result.minimizer, LBFGS())
# optim_result = optimize(planet_model_objective, optim_result.minimizer, LBFGS())

# T0, Pb, e, x, omega0, F0, F1, F2  = optim_result.minimizer
# T0 += T_EP
# Коэффициенты
T0, Pb, e, x, omega0, F0, F1, F2 = fit_result.param
T0 += T_EP
dT0, dPb, de, dx, domega0, dF0, dF1, dF2 = standard_errors(fit_result)
println("Эпоха планеты (T0): $T0 ($dT0)")
println("Период планеты (Pb): $Pb ($Pb)")
println("Эксцентриситет планеты: $(e) ($de)")
println("Проекция большой полуоси планеты (x): $x ($dx)")
println("Угол периастра планеты (omega0): $omega0 ($domega0)")
println("F0: $F0 ($dF0)")
println("F1: $F1 ($dF1)")
println("F1: $F2 ($dF1)")

# Предсказанные значения
predicted_trend = F0 .+ F1 .* t .+ F2 .* t .^2
predicted_signal = planet_model(t, fit_result.param)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
predicted_trend_fine = F0 .+ F1 .* t_fine .+ F2 .* t_fine .^2
predicted_signal_fine = planet_model(t_fine, fit_result.param)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fitted planet")
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

#-----------------------------------------------

detrended_y = y .- predicted_trend
residual_signal = y .- predicted_signal
fig, ax = subplots()
errorbar(t, residual_signal, yerr = y_err, fmt=".")

using LombScargle
fig, ax = subplots()
pgram = lombscargle(t, detrended_y, y_err, maximum_frequency = (N_iters - 1) / (time_finish - time_start))
plot(periodpower(pgram)...)
pgram = lombscargle(t, residual_signal, y_err, maximum_frequency = (N_iters - 1) / (time_finish - time_start))
plot(periodpower(pgram)...)
yscale("log")
xscale("log")

#-----------------------------------------------
param_name = :F1
t = params_values[:PEPOCH] .- T_EP
y = params_values[param_name]
y_err = params_errors[param_name]
y_err_median = median(y_err)
ind = y_err .< Inf * y_err_median

t = t[ind]
y = y[ind]
y_err = y_err[ind]
weights = (y_err) .^ (-2.0)
poly_deg = 1

p = Polynomials.fit(t, y, poly_deg, weights = weights)  # Полином степени 2

# Предсказанный тренд
trend = p.(t)

# Данные без тренда
detrended_y = y .- trend
detrended_y_abs_median = median(abs.(detrended_y))

fig, ax = subplots()
errorbar(t, detrended_y, yerr = y_err, fmt=".")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

#-----------------------------------------------

function v_x(T, T0, Pb, e, x, omega0)
    Pb, e, x = abs.((Pb, e, x))
    M = 2 * pi / Pb * (T - T0)
    U_func = U -> U - e * sin(U) - M
    U = find_zero(U_func, M)
    omega = omega0 * pi / 180.0
    dU_dt = 2 * pi / (Pb * 86400.0) / (1 - e * cos(U)) 
    dRoe = x * sin(omega) * (cos(U) - e) + x * sqrt(1 - e^2) * cos(omega) * sin(U)
    v_x = x * sin(omega) * sin(U) * dU_dt - x * sqrt(1 - e^2) * cos(omega) * cos(U) * dU_dt
    return v_x
end

function a_x(T, T0, Pb, e, x, omega0)
    Pb, e, x = abs.((Pb, e, x))
    M = 2 * pi / Pb * (T - T0)
    U_func = U -> U - e * sin(U) - M
    U = find_zero(U_func, M)
    omega = omega0 * pi / 180.0
    dU_dt = 2 * pi / (Pb * 86400.0) / (1 - e * cos(U)) 
    d2U_dt2 = -dU_dt^2 * e * sin(U) / (1 - e * cos(U))
    dRoe = x * sin(omega) * (cos(U) - e) + x * sqrt(1 - e^2) * cos(omega) * sin(U)
    v_x = x * sin(omega) * sin(U) * dU_dt - x * sqrt(1 - e^2) * cos(omega) * cos(U) * dU_dt
    a_x = x * sin(omega) * (cos(U) * dU_dt^2 + sin(U) * d2U_dt2) + x * sqrt(1 - e^2) * cos(omega) * (sin(U) * dU_dt^2 - cos(U) * d2U_dt2)
    return a_x
end

using LsqFit
using Printf
# Нелинейная модель
planet_model(t, p) = (1.0 .+ v_x.(t, p[1], p[2], p[3], p[4], p[5])) .* (p[6] .+ p[7] .* t) .+ a_x.(t, p[1], p[2], p[3], p[4], p[5]) * F0

# Начальные параметры: [A, omega, T0, offset]
p0 = [T0 - T_EP, Pb, e, x, omega0, p[0], p[1]]

# Подгонка
fit_result = curve_fit(planet_model, t, y, weights, p0)

# Коэффициенты
T0, Pb, e, x, omega0, F1, F2  = fit_result.param
T0 += T_EP
dT0, dPb, de, dx, domega0, dF1, dF2 = standard_errors(fit_result)
println("Эпоха планеты (T0): $T0 ($dT0)")
println("Период планеты (Pb): $Pb ($Pb)")
println("Эксцентриситет планеты: $(e) ($de)")
println("Проекция большой полуоси планеты (x): $x ($dx)")
println("Угол периастра планеты (omega0): $omega0 ($domega0)")
println("F1: $F1 ($dF1)")
println("F1: $F2 ($dF1)")

predicted_trend = p[0] .+ p[1] .* t
predicted_signal = planet_model(t, p0)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
predicted_trend_fine = p[0] .+ p[1] .* t_fine
predicted_signal_fine = planet_model(t_fine, p0)


# Предсказанные значения
predicted_trend = F1 .+ F2 .* t
predicted_signal = planet_model(t, fit_result.param)

t_fine = collect(LinRange(time_start, time_finish, 1000)) .- T_EP
predicted_trend_fine = F1 .+ F2 .* t_fine
predicted_signal_fine = planet_model(t_fine, fit_result.param)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
errorbar(t .+ T_EP, y .- predicted_trend, yerr = y_err, fmt=".")
plot(t_fine .+ T_EP, predicted_signal_fine .- predicted_trend_fine, label="Fitted planet")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

seconds_per_day = 86400.0

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("PL_T0  = %.4f ± %.4f", T0, dT0), "\n",
    @sprintf("PL_PB  = %.2f ± %.2f", Pb, dPb), "\n",
    @sprintf("PL_ECC = %.5f ± %.5f",  e, de), "\n",
    @sprintf("PL_X   = %.5f ± %.5f",  x, dx), "\n",
    @sprintf("PL_OM  = %.2f ± %.2f",  omega0, domega0), "\n",
    @sprintf("F1     = %.8e ± %.2e", F1, dF1), "\n",
    @sprintf("F2     = %.5e ± %.2e", F2 / seconds_per_day, dF2 / seconds_per_day)
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