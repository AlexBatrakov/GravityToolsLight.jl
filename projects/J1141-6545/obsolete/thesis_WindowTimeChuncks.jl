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

#-------------------------------------------------------------------------------------

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551

N_window = 1
N_iters = N_window * 64 + 1
time_arr = collect(LinRange(time_start, time_finish, N_iters))

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_TC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_TC.par",
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

# function extract_parameters(results_global_iters, params_names)
#     N = length(results_global_iters.all_global_iterations)
#     params_values = Dict{Symbol, Vector{Float64}}()
#     params_errors = Dict{Symbol, Vector{Float64}}()
#     for param_name in params_names
#         params_values[param_name] = Float64[]
#         params_errors[param_name] = Float64[]
#     end
#     for i in 1:N
#         try 
#             for param_name in params_names
#                 param_ind = results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters_order[param_name]
#                 push!(params_values[param_name], results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].post_fit)
#                 push!(params_errors[param_name], results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].uncertainty)
#             end
#         catch error
#                 # param_values[i] = NaN
#                 # param_errors[i] = NaN
#         end
#     end
#     return params_values, params_errors
# end

params_values = Dict{Symbol, Vector{Float64}}()
params_errors = Dict{Symbol, Vector{Float64}}()

# param_names = [:PEPOCH, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :M2, :MTOT, :I, :offset, :AR, :OMDOT, :PBDOT, :GAMMA]
# param_names = [:PEPOCH, :F0, :F1, :F2, :offset]
param_names = [:PEPOCH, :F0, :F1, :PB, :T0, :A1, :OM, :ECC, :I, :offset]

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

for param_name in param_names[2:end]
    epoch = 54000.0
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
T_EP = 54000.0
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

P1_estimate = 3450.0
P2_estimate = 1970.0

#-----------------------------------------------

using LsqFit
using Printf
# Нелинейная модель
one_harminic_model(t, p) = p[1] .* cos.(p[2] .* (t .- p[3])) .+ p[4] .+ p[5] .* t .+ p[6] .* t .^ 2

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
legend()
title(L"\Delta F = F - (F_0 + F_1 (T - T_0) + \frac{1}{2} F_2 (T - T_0)^2)")
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


#-----------------------------------------------

# # Нелинейная модель
# sinusoid_model(t, p) = p[1] .* sin.(p[2] .* (t .- p[3])) + p[4] .* sin.(p[5] .* (t .- p[6])) .+ p[7]

# # Начальные параметры: [A1, omega1, phi1, A2, omega2, phi2, offset]
# p0 = [detrended_y_abs_median, 2*pi/P1_estimate, 55500.0, detrended_y_abs_median, 2*pi/P2_estimate, 55500, 0.0]

# # Подгонка
# fit_result = curve_fit(sinusoid_model, t, detrended_y, weights, p0)

# # Коэффициенты
# A1, omega1, T01, A2, omega2, T02, offset = fit_result.param
# dA1, domega1, dT01, dA2, domega2, dT02, doffset = standard_errors(fit_result)
# P1 = 2*pi/omega1
# dP1 = abs(- 2*pi/omega1^2 * domega1)
# P2 = 2*pi/omega2
# dP2 = abs(- 2*pi/omega2^2 * domega2)
# println("1. Амплитуда косинуса (A1): $A1 ($dA1)")
# println("1. Угловая частота (omega1): $omega1 ($domega1)")
# println("1. Период (P1): $(P1) ($dP1)")
# println("1. Время нулевой фазы (T01): $T01 ($dT01)")
# println("2. Амплитуда косинуса (A2): $A2 ($dA2)")
# println("2. Угловая частота (omega2): $omega2 ($domega2)")
# println("2. Период (P2): $(P2) ($dP2)")
# println("2. Время нулевой фазы (T02): $T02 ($dT02)")
# println("Оффсет (offset): $offset ($doffset)")

# # Предсказанные значения
# predicted_sinusoid = sinusoid_model(t, fit_result.param)
# predicted_sinusoid_fine = sinusoid_model(t_fine, fit_result.param)

# fig, ax = subplots()
# errorbar(t, detrended_y, yerr = y_err, fmt=".")
# plot(t_fine, predicted_sinusoid_fine, label="Подогнанные 2 синусоиды")
# ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

# # Форматируем строку с учётом ошибок
# params_text = string(
#     @sprintf("A1 = %.2e ± %.2e", A1, dA1), "\n",
#     @sprintf("P1 = %.2f ± %.2f", P1, dP1), "\n",
#     @sprintf("T01 = %.2f ± %.2f", T01, dT01), "\n",
#     @sprintf("A2 = %.2e ± %.2e", A2, dA2), "\n",
#     @sprintf("P2 = %.2f ± %.2f", P2, dP2), "\n",
#     @sprintf("T02 = %.2f ± %.2f", T02, dT02)
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

# # Оформление
# xlabel("Time (t)")
# ylabel("Signal (y)")
# legend()
# title("Подгонка 2 синусоид с LsqFit")
# grid()
# tight_layout()
# show()

# #-----------------------------------------------

# residual_signal = detrended_y .- predicted_sinusoid
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
t = params_values[:PEPOCH]
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


#-----------------------------------------------




