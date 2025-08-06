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

using Distributed
addprocs(8)

@everywhere using Revise
@everywhere using GravityToolsLight

#-------------------------------------------------------------------------------------

backend_name = "Medusa"

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/backends",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_$(backend_name).par",
    tim_file = "$(backend_name).tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("NITS", 1), TP("TNRedC", 16)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

tim_file_path = joinpath(settings.work_dir, settings.tim_file)
backends = settings.backends

# Чтение данных из файла .tim
tim_file_data = readdlm(tim_file_path, String)
times = parse.(Float64, tim_file_data[3:end, 3])

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = times[1]
time_finish = times[end]
time_validation = time_start + 0.8 * (time_finish - time_start)
# TNRedFlow = log10((time_validation - time_start) / (time_finish - time_start))


global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 6,
    nits = [10,5,5,5,5,5],
    gain = [1,1,1,1,1,1],
    tparams_local = [
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_validation, flag=1)],
        [TP("FINISH", time_finish)]
        ], 
    flags = ["", "", "", "", "", "-nofit"],
    fit_EFACs_EQUADs = [true, true, true, true, true, false],
    print_output = true
)


# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -12.0, max = -8.0, N = 8, range_rule=:lin),
    Var(name = "TNRedGam", min = 1.0, max = 7.0, N = 8, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, max = 1000.0, from_min=true),
    DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_3sigma], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)


run_tempo_general(tf)

for i in 1:3

run_tempo_general(tf, just_refine=true)

chisqr_min = tf.grid.min[:chisqr]

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), rasterized=true)
cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
# plot([], [], label=L"\Delta\chi^{2} = 1\sigma", "-",  color="black")
# plot([], [], label=L"\Delta\chi^{2} = 2\sigma", "--", color="black")
# plot([], [], label=L"\Delta\chi^{2} = 3\sigma", "-.", color="black")
ax.set_ylabel(tf.test_params.x.name, size=14)
ax.set_xlabel(tf.test_params.y.name, size=14)
#title(L"$\chi^2_\mathrm{min} = %$chisqr_min$")
# Добавление плашки с текстом
annot_text = L"$\chi^2_\mathrm{min} = %$chisqr_min$"
# Добавление текста с LaTeX
text(
    0.05, 0.05, annot_text, transform=gca().transAxes,
    fontsize=12,
    bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$")
#ax.invert_xaxis()
# xlim((4.1,5.1))
# ylim((-10.9,-10.5))
# xlim((1.5,2.5))
# ylim((-10.8,-10.2))
title("$(backend_name)")
legend(fontsize=12)
tight_layout()

end

# save("saves/tf_GR_Medusa_20percent.jld", "tf", tf)
# tf = load("saves/tf_GR_Medusa_20percent.jld", "tf")

# tf = load("saves/tf_GR_PL_RN_100d.jld", "tf")

#-------------------------------------------------------------------------------------