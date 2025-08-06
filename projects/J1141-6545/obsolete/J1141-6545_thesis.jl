using PyPlot
using Revise
using DelimitedFiles
using GravityToolsLight
using MultivariateStats
using Statistics
using GLM
using DataFrames
using JLD
pygui(true)

using Distributed
addprocs(8)

@everywhere using Revise
@everywhere using GravityToolsLight

#-------------------------------------------------------------------------------------

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_RN.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 1), TP("TNRedC", 30)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
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
    iters = 2,
    nits = [8, 2],
    gain = [1, 1],
    tparams_local = [
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_finish)]
        ], 
    flags = ["", "-nofit"],
    fit_EFACs_EQUADs = [false, false],
    print_output = true
)


# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -11.0, max = -10.0, N = 7, range_rule=:lin),
    Var(name = "TNRedGam", min = 1.0, max = 5.0, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT),
    desired_refinement_level = 2,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true)
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    # DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_3sigma], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

# run_tempo_general(tf, just_refine=true)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
# plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="black")
# plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="black")
# plot([], [], label=L"\Delta\chi^{2} = 3\sigma", "-.", color="black")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
# xlim((4.1,5.1))
# ylim((-10.9,-10.5))
# xlim((1.5,2.5))
# ylim((-10.8,-10.2))
tight_layout()

# save("saves/tf_GR_RN_1y.jld", "tf", tf)
# tf = load("saves/tf_GR_RN_1y.jld", "tf")

# tf = load("saves/tf_GR_PL_RN_100d.jld", "tf")

#-------------------------------------------------------------------------------------

TNRedC_arr = collect(10:1:120)
N_iters = length(TNRedC_arr)

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_PL_RN_testC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PL_RN.par",
    tim_file = "../J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 5), TP("TNRedAmp", -10.5), TP("TNRedGam", 2.0)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = [i < 25 ? 10 : 5 for i in TNRedC_arr],
    gain = [1],
    tparams_local = [[TP("TNRedC", TNRedC)] for TNRedC in TNRedC_arr],
    # flags = [""],
    # fit_EFACs_EQUADs = [false],
    print_output = true
)

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

#-------------------------------------------------------------------------------------

TNRedC_arr = collect(10:1:120)
N_iters = length(TNRedC_arr)

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_RN_testC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_RN.par",
    tim_file = "../J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 5)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = [5],
    gain = [1],
    tparams_local = [[TP("TNRedC", TNRedC)] for TNRedC in TNRedC_arr],
    # flags = [""],
    # fit_EFACs_EQUADs = [false],
    print_output = true
)

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

#-------------------------------------------------------------------------------------

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551

N_iters = 24
time_arr = collect(LinRange(time_start, time_finish, N_iters))

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_TC",
    version = Tempo2(),
    par_file_init = "../DDSTG_GR_TC.par",
    tim_file = "../J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals -epoch center",
    tparams = [TP("NITS", 20), TP("GAIN_ONE_AFTER_NITS", 15)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = N_iters - 1,
    nits = [20],
    gain = [0.2],
    tparams_local = [[TP("START", time_arr[i], flag=1), TP("FINISH", time_arr[i+1], flag=1)] for i in 1:N_iters-1],
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
            if param_name == :offset
                push!(param_values, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.offset[1])
                push!(param_errors, results_global_iters.all_global_iterations[i].last_internal_iteration.result.basic.offset[2])
            else
                param_ind = results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters_order[param_name]
                push!(param_values, results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].post_fit)
                push!(param_errors, results_global_iters.all_global_iterations[i].last_internal_iteration.result.fit_parameters[param_ind].uncertainty)
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

for param_name in [:PEPOCH, :F0, :F1, :PB, :A1, :OM, :ECC, :I, :offset]
    params_values[param_name], params_errors[param_name] = extract_parameter(results_global_iters, param_name)
end

# params_values, params_errors = extract_parameters(results_global_iters, [:PEPOCH, :F0, :F1, :PB, :A1, :OM, :ECC, :I])

# PEPOCH_values, _ = extract_parameter(results_global_iters, :PEPOCH)
# F0_values, F0_errors = extract_parameter(results_global_iters, :F0)
# F1_values, F1_errors = extract_parameter(results_global_iters, :F1)
# PB_values, PB_errors = extract_parameter(results_global_iters, :PB)
# A1_values, A1_errors = extract_parameter(results_global_iters, :A1)
# OM_values, OM_errors = extract_parameter(results_global_iters, :OM)
# ECC_values, ECC_errors = extract_parameter(results_global_iters, :ECC)

fig, ax = subplots()
errorbar(params_values[:PEPOCH], params_values[:F0], params_errors[:F0])



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

using Polynomials
param_name = :F1
t = params_values[:PEPOCH]
y = params_values[param_name]
y_err = params_errors[param_name]
y_err_median = median(y_err)
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

#-----------------------------------------------

using LsqFit
using Printf
# Нелинейная модель
sinusoid_model(t, p) = p[1] .* sin.(p[2] .* (t .- p[3])) .+ p[4]

# Начальные параметры: [A, omega, T0, offset]
p0 = [detrended_y_abs_median, 2*pi/3400.0, 54400.0, 0.0]

# Подгонка
fit_result = curve_fit(sinusoid_model, t, detrended_y, weights, p0)

# Коэффициенты
A, omega, T0, offset = fit_result.param
dA, domega, dT0, doffset = standard_errors(fit_result)
P = 2*pi/omega
dP = abs(- 2*pi/omega^2 * domega)
println("Амплитуда косинуса (A): $A ($dA)")
println("Угловая частота (omega): $omega ($domega)")
println("Период (P): $(P) ($dP)")
println("Время нулевой фазы (T0): $T0 ($dT0)")
println("Оффсет (offset): $offset ($doffset)")

# Предсказанные значения
predicted_simusoid = sinusoid_model(t, fit_result.param)

fig, ax = subplots()
errorbar(t, detrended_y, yerr = y_err, fmt=".")
plot(t, predicted_simusoid, label="Подогнанная синусоида")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("A = %.2e ± %.2e", A, dA), "\n",
    @sprintf("P = %.2f ± %.2f", P, dP), "\n",
    @sprintf("T0 = %.2f ± %.2f", T0, dT0)
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

# Нелинейная модель
sinusoid_model(t, p) = p[1] .* sin.(p[2] .* (t .- p[3])) + p[4] .* sin.(p[5] .* (t .- p[6])) .+ p[7]

# Начальные параметры: [A1, omega1, phi1, A2, omega2, phi2, offset]
p0 = [detrended_y_abs_median, 2*pi/3400.0, 54500.0, detrended_y_abs_median, 2*pi/2100.0, 54500, 0.0]

# Подгонка
fit_result = curve_fit(sinusoid_model, t, detrended_y, weights, p0)

# Коэффициенты
A1, omega1, T01, A2, omega2, T02, offset = fit_result.param
dA1, domega1, dT01, dA2, domega2, dT02, doffset = standard_errors(fit_result)
P1 = 2*pi/omega1
dP1 = abs(- 2*pi/omega1^2 * domega1)
P2 = 2*pi/omega2
dP2 = abs(- 2*pi/omega2^2 * domega2)
println("1. Амплитуда косинуса (A1): $A1 ($dA1)")
println("1. Угловая частота (omega1): $omega1 ($domega1)")
println("1. Период (P1): $(P1) ($dP1)")
println("1. Время нулевой фазы (T01): $T01 ($dT01)")
println("2. Амплитуда косинуса (A2): $A2 ($dA2)")
println("2. Угловая частота (omega2): $omega2 ($domega2)")
println("2. Период (P2): $(P2) ($dP2)")
println("2. Время нулевой фазы (T02): $T02 ($dT02)")
println("Оффсет (offset): $offset ($doffset)")

# Предсказанные значения
predicted_simusoid = sinusoid_model(t, fit_result.param)

fig, ax = subplots()
errorbar(t, detrended_y, yerr = y_err, fmt=".")
plot(t, predicted_simusoid, label="Подогнанные 2 синусоиды")
ylim((-5*detrended_y_abs_median, 5*detrended_y_abs_median))

# Форматируем строку с учётом ошибок
params_text = string(
    @sprintf("A1 = %.2e ± %.2e", A1, dA1), "\n",
    @sprintf("P1 = %.2f ± %.2f", P1, dP1), "\n",
    @sprintf("T01 = %.2f ± %.2f", T01, dT01), "\n",
    @sprintf("A2 = %.2e ± %.2e", A2, dA2), "\n",
    @sprintf("P2 = %.2f ± %.2f", P2, dP2), "\n",
    @sprintf("T02 = %.2f ± %.2f", T02, dT02)
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
title("Подгонка 2 синусоид с LsqFit")
grid()
tight_layout()
show()