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

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_fake.par",
    tim_file = "DDSTG_GR_fake.simulate",
    # par_file_init = "DDSTG_GR_ABE_test.par",
    # tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("STG_ALPHA0", 0.0), TP("STG_BETA0", 0.0), TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )
# bsets = basic_settings

# results_basic = run_tempo_basic(basic_settings)

# extract_internal_iterations_values(results_basic, :MTOT, :post_fit)

# global_iters_settings = GlobalIterationsSettings(
#     keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
#     iters = 2,
#     nits = [2,2],
#     gain = [1.0, 1.0],
#     tparams_local = [
#         [TP("P_LAMBDA", flag=0), TP("P_ETA", flag=0)],
#         [TP("P_LAMBDA", flag=1), TP("P_ETA", flag=1)]
#         ]
#     )
# gisets = global_iters_settings

# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)




# general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)


general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "P_PHI",    min = 0.0, max = 360.0, N = 7, range_rule=:lin),
    Var(name = "P_DELTA", min = 0.0, max = 180.0, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :EDOT, :PBDOT, :XDOT, :X2DOT, :OMDOT, :GAMMA, :I, :M2, :MTOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_LAMBDADOT, :P_LAMBDA2DOT, :P_ETA, :P_ETADOT, :P_ETA2DOT, :P_DELTA, :P_PHI, :P_PHIDOT, :ABE_X, :ABE_XDOT, :ABE_X2DOT, :ABE_M, :ABE_PB, :ABE_PBDOT, :ABE_E, :ABE_EDOT),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true),
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:chisqr, diffs = [1.0, 1.0], contours = [lvl_1sigma, lvl_2sigma], from_min=true)
    )


tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

run_tempo_general(tf, just_refine=true)

chisqr_min = minimum(tf.grid.vars[:chisqr])

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="red")
plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="red")
plot([], [], label=L"\Delta\chi^{2} = 3\sigma", "-.", color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

local_mins = Tuple{Int64, Int64}[]
local_mins_mask = BitMatrix(undef, tf.grid.x.N, tf.grid.y.N)

cell_selector_status = fill(false, tf.grid.x.N, tf.grid.y.N)
for i_cell in 1:tf.grid.x.N, j_cell in 1:tf.grid.y.N
    if GravityToolsLight.cell_selector(i_cell, j_cell, tf.grid, tf.ref_sets.units[1], at_corner=true)
        push!(local_mins, (i_cell, j_cell))
        local_mins_mask[i_cell,j_cell] = 1
        # println("$(tf.grid.y.values[j_cell]), $(tf.grid.x.values[i_cell]), $(tf.grid.vars[:chisqr][i_cell,j_cell])")
        cell_selector_status[i_cell, j_cell] = true
        plot(tf.grid.y.values[j_cell], tf.grid.x.values[i_cell], ".", color="black")
        plt.annotate("$(tf.grid.vars[:chisqr][i_cell,j_cell])", (tf.grid.y.values[j_cell], tf.grid.x.values[i_cell]), textcoords="offset points", xytext=(1,1), ha="center")
    end
end

plot(78.080553383129876067, 3.3629878534675309042, "x", color="violet")

cells_to_refine = 0
cell_selector_status = fill(false, tf.grid.x.N, tf.grid.y.N)
for i_cell in 1:tf.grid.x.N, j_cell in 1:tf.grid.y.N
    if GravityToolsLight.cell_selector(i_cell, j_cell, tf.grid)
        # println("$(tf.grid.y.values[j_cell]), $(tf.grid.x.values[i_cell]), $(tf.grid.vars[:chisqr][i_cell,j_cell])")
        cell_selector_status[i_cell, j_cell] = true
        cells_to_refine +=1
        plot(tf.grid.y.values[j_cell], tf.grid.x.values[i_cell], ".", color="black")
        plot(tf.grid.y.values[j_cell+1], tf.grid.x.values[i_cell], ".", color="black")
        plot(tf.grid.y.values[j_cell], tf.grid.x.values[i_cell+1], ".", color="black")
        plot(tf.grid.y.values[j_cell+1], tf.grid.x.values[i_cell+1], ".", color="black")
        # plt.annotate("$(tf.grid.vars[:chisqr][i_cell,j_cell])", (tf.grid.y.values[j_cell], tf.grid.x.values[i_cell]), textcoords="offset points", xytext=(1,1), ha="center")
    end
end

open("local_mins.dat", "w") do file
    # for (i, j) in local_mins
        # write(file, "P_PHI = $(tf.grid.y.values[j]), P_DELTA = $(tf.grid.x.values[i])\n")
        ind_global_min = findmin(tf.grid.vars[:dchisqr])[2]
        for key in keys(tf.grid.vars)
            write(file, "delta $key = $(round.(tf.grid.vars[key][local_mins_mask] .- tf.grid.vars[key][ind_global_min], sigdigits=5))\n")
        end
        write(file, "\n\n\n")

        for key in keys(tf.grid.vars)
            write(file, "$key = $((tf.grid.vars[key][local_mins_mask]))\n")
        end
    # end
end

for i_cell in 1:tf.grid.x.N, j_cell in 1:tf.grid.y.N
    if tf.grid.status[i_cell, j_cell] == 1
        plot(tf.grid.y.values[j_cell], tf.grid.x.values[i_cell], ".", color="black")
    end
end

var = :X2DOT
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[var], cmap="Blues_r", rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="red")
plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="red")
plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "-.", color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

fig, ax = subplots()
mask = tf.grid.status .== 1 .&& (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr] .< 50)
plot(tf.grid.vars[:P_LAMBDA][mask], tf.grid.vars[:chisqr][mask] .- tf.grid.min[:chisqr], ".", color="black")
ylim((0, 30))

tf.grid.vars[:dchisqr] = tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]

fig, ax = subplots()
scatter(tf.grid.vars[:AFACDOT][mask], tf.grid.vars[:BFACDOT][mask], c = tf.grid.vars[:chisqr][mask] .- tf.grid.min[:chisqr] )
colorbar()

tf.grid.vars[:P_PHI]   = [P_PHI for P_PHI in tf.grid.x.values, P_DELTA in tf.grid.y.values]
tf.grid.vars[:P_DELTA] = [P_DELTA for P_PHI in tf.grid.x.values, P_DELTA in tf.grid.y.values]

# Список ключей, который ты дал
keys = [:dchisqr, :P_PHI, :P_DELTA, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :GAMMA, :I, :M2, :MTOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_ETA]

keys = [:dchisqr, :P_PHI, :P_DELTA, :PB, :T0, :OM, :ECC, :XDOT, :X2DOT, :OMDOT, :I, :M2, :MTOT]

keys = [:dchisqr, :P_ETA, :P_LAMBDA, :PB, :T0, :OM, :ECC, :XDOT, :X2DOT, :OMDOT, :I, :M2, :MTOT]

# Создаем пустую матрицу для данных
data_matrix = []

# Цикл по ключам
for key in keys
    # Линеаризуем данные для каждого ключа
    data = tf.grid.vars[key][mask]  # Это предполагает, что mask уже применен
    linear_data = vec(data)  # Преобразование в одномерный массив
    push!(data_matrix, linear_data)
end

# Преобразуем в матрицу, где строки — это наблюдения, а столбцы — переменные
data_matrix = hcat(data_matrix...)

# Функция для стандартизации данных: вычитание среднего и деление на стандартное отклонение
function standardize_data(data)
    return (data .- mean(data, dims=1)) ./ std(data, dims=1)
end

# Стандартизация матрицы данных
data_standardized = standardize_data(data_matrix)

# Применение PCA
model = fit(PCA, data_standardized', pratio=1)  # Транспонируем, так как PCA ожидает переменные по столбцам

# Получаем матрицу коэффициентов главных компонент
coefficients = model.proj  # Это коэффициенты разложения главных компонент

# Исправляем построение тепловой карты (оси должны быть правильными)
figure()
imshow(coefficients', cmap="RdBu", aspect="auto")  # Транспонируем коэффициенты для правильного отображения осей
colorbar()

# Указываем подписи по осям
xticks(0:length(keys)-1, string.(keys), rotation=45)  # Переменные по горизонтали
yticks(0:size(coefficients, 2)-1, ["PC$(i)" for i in 1:size(coefficients, 2)])  # Главные компоненты по вертикали

xlabel("Variables")
ylabel("Principal Components")
title("Heatmap of PCA Coefficients (Eigenvectors)")
show()

# Теперь строим график объясненной дисперсии каждой PCA

# Получаем собственные значения (дисперсии) для каждой компоненты
explained_variances = model.prinvars

# Нормализуем значения, чтобы выразить их в процентах от общей дисперсии
explained_variance_ratio = explained_variances / sum(explained_variances) * 100

# Строим график
figure()
bar(1:length(explained_variance_ratio), explained_variance_ratio)
xlabel("Principal Component")
ylabel("Explained Variance (%)")
title("Explained Variance by Principal Components")
show()

# Построим корреляционную матрицу для исходных переменных
correlation_matrix = cor(data_standardized)

# Построение тепловой карты для матрицы корреляций
figure()
imshow(correlation_matrix, cmap="RdBu", aspect="auto")
colorbar()

# Указываем подписи по осям
xticks(0:length(keys)-1, string.(keys), rotation=45)  # Переменные по горизонтали
yticks(0:length(keys)-1, string.(keys))  # Переменные по вертикали

xlabel("Variables")
ylabel("Variables")
title("Correlation Matrix of Original Variables")
show()

# save("tf_ABE.jld", "tf", tf)
tf = load("tf_ABE.jld", "tf")

# save("tf_ABE_fIDOT.jld", "tf", tf)
tf = load("tf_ABE_fIDOT.jld", "tf")

# save("tf_ABE_fIDOT_I2DOT.jld", "tf", tf)
tf = load("tf_ABE_fIDOT_I2DOT.jld", "tf")

# save("tf_ABE_fIDOT_I2DOT_OMDOT.jld", "tf", tf)
tf = load("tf_ABE_fIDOT_I2DOT_OMDOT.jld", "tf")
#-------------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    # par_file_init = "DDSTG_GR_fake.par",
    # tim_file = "DDSTG_GR_fake.simulate",
    par_file_init = "DDSTG_GR_PREC_test.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("STG_ALPHA0", 0.0), TP("STG_BETA0", 0.0), TP("NITS", 5)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "C_XI",    min = 0.0, max = 3.0, N = 7, range_rule=:lin),
    Var(name = "C_DELTA", min = 0.0, max = 180.0, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
    )

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :EDOT, :PBDOT, :XDOT, :XDOT_I_S, :XDOT_I_Q, :X2DOT, :GAMMA, :M2, :MTOT, :J_I, :J_DELTA, :C_XI, :C_DELTA, :PREC_OM_SEC, :PREC_OMDOT_SEC, :PREC_OMDOT_PN_SEC, :PREC_PHI, :PREC_PSI),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true),
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:chisqr, diffs = [1.0, 1.0], contours = [lvl_1sigma, lvl_2sigma], from_min=true)
    )


tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

#-------------------------------------------------------------------------------------



basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_test_RN.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 1), TP("TNRedC", 83)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basic
# EFACs, EQUADs, log10EQUADs, obj = calculate_EFACs_EQUADs(basic_settings)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 5,
    nits = [4,4,4,4,4],
    gain = [1.0, 1.0, 1.0, 1.0, 1.0]
    )

# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -13, max = -9, N = 5, range_rule=:lin),
    Var(name = "TNRedGam", min = 1, max = 5, N = 5, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:res_obj, :chisqr, :rms_post_fit, :rms_tn_post_fit, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :GAMMA, :I, :M2, :MTOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_ETA),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true),
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:chisqr, diffs = [1.0], contours = [1.0], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=1000), rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="red")
plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="red")
plot([], [], label=L"\Delta\chi^{2} = 3\sigma", "-.", color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

Amp_arr = [10.0 ^ TNRedAmp for TNRedAmp in tf.grid.x.values, TNRedGam in tf.grid.y.values]
lambda = 1e10

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr] .+ lambda .* Amp_arr, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_3sigma), rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr] .+ lambda .* Amp_arr, levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="red")
plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="red")
plot([], [], label=L"\Delta\chi^{2} = 3\sigma", "-.", color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}_{reg}$")
#ax.invert_xaxis()
tight_layout()

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_tn_post_fit] .- tf.grid.min[:rms_tn_post_fit], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100.0), rasterized=true)
# cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_tn_post_fit] .- tf.grid.min[:rms_tn_post_fit], levels=[1.0], linestyles=["-"], colors="red")
# plot([], [], label=L"\Delta RMS = 1", "-",  color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$RMS_TN - min(RMS_TN)$")
#ax.invert_xaxis()
tight_layout()


rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_tn_post_fit] .- tf.grid.min[:rms_tn_post_fit], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100.0), rasterized=true)
# cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_tn_post_fit] .- tf.grid.min[:rms_tn_post_fit], levels=[1.0], linestyles=["-"], colors="red")
# plot([], [], label=L"\Delta RMS = 1", "-",  color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$RMS_TN - min(RMS_TN)$")
#ax.invert_xaxis()
tight_layout()


rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_post_fit] .- tf.grid.min[:rms_post_fit], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10000), rasterized=true)
# cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:rms_post_fit] .- tf.grid.min[:rms_post_fit], levels=[1.0], linestyles=["-"], colors="red")
# plot([], [], label=L"\Delta RMS = 1", "-",  color="red")
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$RMS - min(RMS)$")
#ax.invert_xaxis()
tight_layout()



#-------------------------------------------------------------------------------------



basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PREC_test.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )


global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 9,
    nits = [5,1,5,1,5,1,5,1,5,1,5],
    gain = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    tparams_local = [
        [TP("C_XI", flag=0), TP("C_DELTA", flag=0)],
        [TP("C_XI", flag=1), TP("C_DELTA", flag=1)],
        [TP("C_XI", flag=0), TP("C_DELTA", flag=0)],
        [TP("C_XI", flag=1), TP("C_DELTA", flag=1)],
        [TP("C_XI", flag=0), TP("C_DELTA", flag=0)],
        [TP("C_XI", flag=1), TP("C_DELTA", flag=1)],
        [TP("C_XI", flag=0), TP("C_DELTA", flag=0)],
        [TP("C_XI", flag=1), TP("C_DELTA", flag=1)],
        [TP("C_XI", flag=0), TP("C_DELTA", flag=0)]
        ]
    )

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

#-------------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_CHUNCK.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )


global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = 23,
    nits = 20 * ones(23),
    gain = 0.1 * ones(23),
    tparams_local = [
        [TP("START", 51630.78234259175, flag=1),  TP("FINISH", 51995.50746543971, flag=1)],
        [TP("START", 51995.50746543971, flag=1),  TP("FINISH", 52360.23258828767, flag=1)],
        [TP("START", 52360.23258828767, flag=1),  TP("FINISH", 52724.95771113563, flag=1)],
        [TP("START", 52724.95771113563, flag=1),  TP("FINISH", 53089.68283398358, flag=1)],
        [TP("START", 53089.68283398358, flag=1),  TP("FINISH", 53454.407956831536, flag=1)],
        [TP("START", 53454.407956831536, flag=1), TP("FINISH", 53819.133079679494, flag=1)],
        [TP("START", 53819.133079679494, flag=1), TP("FINISH", 54183.85820252745, flag=1)],
        [TP("START", 54183.85820252745, flag=1),  TP("FINISH", 54548.58332537541, flag=1)],
        [TP("START", 54548.58332537541, flag=1),  TP("FINISH", 54913.30844822337, flag=1)],
        [TP("START", 54913.30844822337, flag=1),  TP("FINISH", 55278.033571071326, flag=1)],
        [TP("START", 55278.033571071326, flag=1), TP("FINISH", 55642.758693919284, flag=1)],
        [TP("START", 55642.758693919284, flag=1), TP("FINISH", 56007.483816767235, flag=1)],
        [TP("START", 56007.483816767235, flag=1), TP("FINISH", 56372.20893961519, flag=1)],
        [TP("START", 56372.20893961519, flag=1),  TP("FINISH", 56736.93406246315, flag=1)],
        [TP("START", 56736.93406246315, flag=1),  TP("FINISH", 57101.65918531111, flag=1)],
        [TP("START", 57101.65918531111, flag=1),  TP("FINISH", 57466.38430815907, flag=1)],
        [TP("START", 57466.38430815907, flag=1),  TP("FINISH", 57831.109431007026, flag=1)],
        [TP("START", 57831.109431007026, flag=1), TP("FINISH", 58195.834553854984, flag=1)],
        [TP("START", 58195.834553854984, flag=1), TP("FINISH", 58560.55967670294, flag=1)],
        [TP("START", 58560.55967670294, flag=1),  TP("FINISH", 58925.28479955089, flag=1)],
        [TP("START", 58925.28479955089, flag=1),  TP("FINISH", 59290.00992239885, flag=1)],
        [TP("START", 59290.00992239885, flag=1),  TP("FINISH", 59654.73504524681, flag=1)],
        [TP("START", 59654.73504524681, flag=1),  TP("FINISH", 60019.46016809477, flag=1)]
        ]
    )

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)





#-------------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_CHUNCK.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar",
    tparams = [TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )


global_iters_settings = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=false, save_global_iterations=true),
    iters = 92,
    nits = 5 * ones(92),
    gain = 1 * ones(92),
    tparams_local = [
        [TP("START", 51630.78234259175, flag=1),  TP("FINISH", 51721.96362330374, flag=1)],
        [TP("START", 51721.96362330374, flag=1),  TP("FINISH", 51813.14490401573, flag=1)],
        [TP("START", 51813.14490401573, flag=1),  TP("FINISH", 51904.32618472772, flag=1)],
        [TP("START", 51904.32618472772, flag=1),  TP("FINISH", 51995.50746543971, flag=1)],
        [TP("START", 51995.50746543971, flag=1),  TP("FINISH", 52086.688746151696, flag=1)],
        [TP("START", 52086.688746151696, flag=1), TP("FINISH", 52177.87002686369, flag=1)],
        [TP("START", 52177.87002686369, flag=1),  TP("FINISH", 52269.051307575675, flag=1)],
        [TP("START", 52269.051307575675, flag=1), TP("FINISH", 52360.23258828767, flag=1)],
        [TP("START", 52360.23258828767, flag=1),  TP("FINISH", 52451.413868999654, flag=1)],
        [TP("START", 52451.413868999654, flag=1), TP("FINISH", 52542.59514971165, flag=1)],
        [TP("START", 52542.59514971165, flag=1),  TP("FINISH", 52633.77643042363, flag=1)],
        [TP("START", 52633.77643042363, flag=1),  TP("FINISH", 52724.95771113563, flag=1)],
        [TP("START", 52724.95771113563, flag=1),  TP("FINISH", 52816.13899184761, flag=1)],
        [TP("START", 52816.13899184761, flag=1),  TP("FINISH", 52907.320272559606, flag=1)],
        [TP("START", 52907.320272559606, flag=1), TP("FINISH", 52998.50155327159, flag=1)],
        [TP("START", 52998.50155327159, flag=1),  TP("FINISH", 53089.68283398358, flag=1)],
        [TP("START", 53089.68283398358, flag=1),  TP("FINISH", 53180.86411469557, flag=1)],
        [TP("START", 53180.86411469557, flag=1),  TP("FINISH", 53272.04539540756, flag=1)],
        [TP("START", 53272.04539540756, flag=1),  TP("FINISH", 53363.22667611955, flag=1)],
        [TP("START", 53363.22667611955, flag=1),  TP("FINISH", 53454.407956831536, flag=1)],
        [TP("START", 53454.407956831536, flag=1), TP("FINISH", 53545.58923754353, flag=1)],
        [TP("START", 53545.58923754353, flag=1),  TP("FINISH", 53636.770518255515, flag=1)],
        [TP("START", 53636.770518255515, flag=1), TP("FINISH", 53727.95179896751, flag=1)],
        [TP("START", 53727.95179896751, flag=1),  TP("FINISH", 53819.133079679494, flag=1)],
        [TP("START", 53819.133079679494, flag=1), TP("FINISH", 53910.31436039149, flag=1)],
        [TP("START", 53910.31436039149, flag=1),  TP("FINISH", 54001.49564110347, flag=1)],
        [TP("START", 54001.49564110347, flag=1),  TP("FINISH", 54092.676921815466, flag=1)],
        [TP("START", 54092.676921815466, flag=1), TP("FINISH", 54183.85820252745, flag=1)],
        [TP("START", 54183.85820252745, flag=1),  TP("FINISH", 54275.039483239445, flag=1)],
        [TP("START", 54275.039483239445, flag=1), TP("FINISH", 54366.22076395143, flag=1)],
        [TP("START", 54366.22076395143, flag=1),  TP("FINISH", 54457.40204466342, flag=1)],
        [TP("START", 54457.40204466342, flag=1),  TP("FINISH", 54548.58332537541, flag=1)],
        [TP("START", 54548.58332537541, flag=1),  TP("FINISH", 54639.764606087396, flag=1)],
        [TP("START", 54639.764606087396, flag=1), TP("FINISH", 54730.94588679939, flag=1)],
        [TP("START", 54730.94588679939, flag=1),  TP("FINISH", 54822.127167511375, flag=1)],
        [TP("START", 54822.127167511375, flag=1), TP("FINISH", 54913.30844822337, flag=1)],
        [TP("START", 54913.30844822337, flag=1),  TP("FINISH", 55004.489728935354, flag=1)],
        [TP("START", 55004.489728935354, flag=1), TP("FINISH", 55095.67100964735, flag=1)],
        [TP("START", 55095.67100964735, flag=1),  TP("FINISH", 55186.85229035933, flag=1)],
        [TP("START", 55186.85229035933, flag=1),  TP("FINISH", 55278.033571071326, flag=1)],
        [TP("START", 55278.033571071326, flag=1), TP("FINISH", 55369.21485178331, flag=1)],
        [TP("START", 55369.21485178331, flag=1),  TP("FINISH", 55460.396132495305, flag=1)],
        [TP("START", 55460.396132495305, flag=1), TP("FINISH", 55551.57741320729, flag=1)],
        [TP("START", 55551.57741320729, flag=1),  TP("FINISH", 55642.758693919284, flag=1)],
        [TP("START", 55642.758693919284, flag=1), TP("FINISH", 55733.93997463127, flag=1)],
        [TP("START", 55733.93997463127, flag=1),  TP("FINISH", 55825.121255343256, flag=1)],
        [TP("START", 55825.121255343256, flag=1), TP("FINISH", 55916.30253605525, flag=1)],
        [TP("START", 55916.30253605525, flag=1),  TP("FINISH", 56007.483816767235, flag=1)],
        [TP("START", 56007.483816767235, flag=1), TP("FINISH", 56098.66509747923, flag=1)],
        [TP("START", 56098.66509747923, flag=1),  TP("FINISH", 56189.846378191214, flag=1)],
        [TP("START", 56189.846378191214, flag=1), TP("FINISH", 56281.02765890321, flag=1)],
        [TP("START", 56281.02765890321, flag=1),  TP("FINISH", 56372.20893961519, flag=1)],
        [TP("START", 56372.20893961519, flag=1),  TP("FINISH", 56463.39022032719, flag=1)],
        [TP("START", 56463.39022032719, flag=1),  TP("FINISH", 56554.57150103917, flag=1)],
        [TP("START", 56554.57150103917, flag=1),  TP("FINISH", 56645.752781751165, flag=1)],
        [TP("START", 56645.752781751165, flag=1), TP("FINISH", 56736.93406246315, flag=1)],
        [TP("START", 56736.93406246315, flag=1),  TP("FINISH", 56828.115343175145, flag=1)],
        [TP("START", 56828.115343175145, flag=1), TP("FINISH", 56919.29662388713, flag=1)],
        [TP("START", 56919.29662388713, flag=1),  TP("FINISH", 57010.47790459912, flag=1)],
        [TP("START", 57010.47790459912, flag=1),  TP("FINISH", 57101.65918531111, flag=1)],
        [TP("START", 57101.65918531111, flag=1),  TP("FINISH", 57192.8404660231, flag=1)],
        [TP("START", 57192.8404660231, flag=1),   TP("FINISH", 57284.02174673509, flag=1)],
        [TP("START", 57284.02174673509, flag=1),  TP("FINISH", 57375.203027447074, flag=1)],
        [TP("START", 57375.203027447074, flag=1), TP("FINISH", 57466.38430815907, flag=1)],
        [TP("START", 57466.38430815907, flag=1),  TP("FINISH", 57557.56558887105, flag=1)],
        [TP("START", 57557.56558887105, flag=1),  TP("FINISH", 57648.74686958305, flag=1)],
        [TP("START", 57648.74686958305, flag=1),  TP("FINISH", 57739.92815029503, flag=1)],
        [TP("START", 57739.92815029503, flag=1),  TP("FINISH", 57831.109431007026, flag=1)],
        [TP("START", 57831.109431007026, flag=1), TP("FINISH", 57922.29071171901, flag=1)],
        [TP("START", 57922.29071171901, flag=1),  TP("FINISH", 58013.471992431005, flag=1)],
        [TP("START", 58013.471992431005, flag=1), TP("FINISH", 58104.65327314299, flag=1)],
        [TP("START", 58104.65327314299, flag=1),  TP("FINISH", 58195.834553854984, flag=1)],
        [TP("START", 58195.834553854984, flag=1), TP("FINISH", 58287.01583456697, flag=1)],
        [TP("START", 58287.01583456697, flag=1),  TP("FINISH", 58378.19711527896, flag=1)],
        [TP("START", 58378.19711527896, flag=1),  TP("FINISH", 58469.37839599095, flag=1)],
        [TP("START", 58469.37839599095, flag=1),  TP("FINISH", 58560.55967670294, flag=1)],
        [TP("START", 58560.55967670294, flag=1),  TP("FINISH", 58651.74095741493, flag=1)],
        [TP("START", 58651.74095741493, flag=1),  TP("FINISH", 58742.922238126914, flag=1)],
        [TP("START", 58742.922238126914, flag=1), TP("FINISH", 58834.10351883891, flag=1)],
        [TP("START", 58834.10351883891, flag=1),  TP("FINISH", 58925.28479955089, flag=1)],
        [TP("START", 58925.28479955089, flag=1),  TP("FINISH", 59016.466080262886, flag=1)],
        [TP("START", 59016.466080262886, flag=1), TP("FINISH", 59107.64736097487, flag=1)],
        [TP("START", 59107.64736097487, flag=1),  TP("FINISH", 59198.828641686865, flag=1)],
        [TP("START", 59198.828641686865, flag=1), TP("FINISH", 59290.00992239885, flag=1)],
        [TP("START", 59290.00992239885, flag=1),  TP("FINISH", 59381.191203110844, flag=1)],
        [TP("START", 59381.191203110844, flag=1), TP("FINISH", 59472.37248382283, flag=1)],
        [TP("START", 59472.37248382283, flag=1),  TP("FINISH", 59563.55376453482, flag=1)],
        [TP("START", 59563.55376453482, flag=1),  TP("FINISH", 59654.73504524681, flag=1)],
        [TP("START", 59654.73504524681, flag=1),  TP("FINISH", 59745.9163259588, flag=1)],
        [TP("START", 59745.9163259588, flag=1),   TP("FINISH", 59837.09760667079, flag=1)],
        [TP("START", 59837.09760667079, flag=1),  TP("FINISH", 59928.27888738278, flag=1)],
        [TP("START", 59928.27888738278, flag=1),  TP("FINISH", 60019.46016809477, flag=1)]
        ]
    )

results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

#-------------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_test_RN.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 4)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basic
# EFACs, EQUADs, log10EQUADs, obj = calculate_EFACs_EQUADs(basic_settings); EFACs

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 5,
    nits = 3 * ones(5),
    gain = 1 * ones(5)
    )

# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -14, max = -8, N = 7, range_rule=:lin),
    Var(name = "TNRedGam", min = 0, max = 6, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:res_obj, :chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :GAMMA, :I, :M2, :MTOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_ETA),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true),
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    DiffContourUnit(:chisqr, diffs = [1.0], contours = [1.0], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10000), rasterized=true)
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

#-------------------------------------------------------------------------------------

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_test_RN.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 1)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basic
# EFACs, EQUADs, log10EQUADs, obj = calculate_EFACs_EQUADs(basic_settings); EFACs

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 3,
    nits = 3 * ones(3),
    gain = 1 * ones(3)
)

# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

# general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

time_chuncks_settings = TimeChuncksSettings(
    number_of_chuncks = 10
)

results_time_chuncks = run_tempo_global_iters(basic_settings, global_iters_settings, time_chuncks_settings)

#-------------------------------------------------------------------------------------

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_final",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_RN.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000  -newpar -writeres -residuals",
    tparams = [TP("NITS", 1)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 2,
    nits = [6, 1],
    gain = [1, 1],
    tparams_local = [
        [TP("TNRedFLow", -0.01933356497235811), TP("FINISH", 59654.21016809477, flag=1)],
        [TP("TNRedFLow", 0.0),                  TP("FINISH", 60019.460168094764551)]
        ], 
    flags = ["", "-nofit"],
    fit_EFACs_EQUADs = [false, false],
    print_output = true
)


# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -13, max = -9, N = 7, range_rule=:lin),
    Var(name = "TNRedGam", min = 1, max = 6, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :GAMMA, :I, :M2, :MTOT, :AFAC, :BFAC, :AFACDOT, :BFACDOT, :AFAC2DOT, :BFAC2DOT, :P_LAMBDA, :P_ETA),
    desired_refinement_level = 4,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true)
    # RelDiffUnit(:chisqr, rel_diff = 1.0, max = 20.0, from_min=true)
#    ContourUnit(:val1, contours = [0.5])
    # DiffContourUnit(:chisqr, diffs = [1.0], contours = [1.0], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

run_tempo_general(tf)

# run_tempo_general(tf, just_refine=true)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), rasterized=true)
ax.set_ylabel(tf.test_params.x.name, size=16)
ax.set_xlabel(tf.test_params.y.name, size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()
