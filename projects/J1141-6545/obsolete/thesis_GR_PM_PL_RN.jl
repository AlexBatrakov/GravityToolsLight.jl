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

#----------------------------------------------------------------------------------
TNRedC = 84

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_GR_PM_PL_RN",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_PM_PL_RN$TNRedC.par",
    tim_file = "J1141-6545_pn_new.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("TNRedC", TNRedC)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
)

# results_basic = run_tempo_basic(basic_settings)
# results_basic.last_internal_iteration.result.basics
# EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings, time_start = 51630.782342591750002, time_finish = 59654.21016809477); EFACs

time_start  = 51630.782342591750002
time_finish = 60019.460168094764551
time_validation = time_finish - 200.0
# TNRedFlow = log10((time_validation - time_start) / (time_finish - time_start))

N_iters = 7

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = N_iters,
    nits = [5,5,5,5,5,5,1],
    gain = ones(N_iters),
    tparams_local = [
        if i < N_iters
            [TP("FINISH", time_validation, flag=1)]
        else
            [TP("FINISH", time_finish)]
        end
        for i in 1:N_iters], 
    flags = [i < N_iters ? "" : "-nofit" for i in 1:N_iters],
    fit_EFACs_EQUADs = [1 <= i < N_iters - 1 ? true : false for i in 1:N_iters],
    print_output = true
)


# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings, global_iters_settings)

test_params = TestParameters(
    Var(name = "TNRedAmp", min = -11.5, max = -9.5, N = 5, range_rule=:lin),
    Var(name = "TNRedGam", min =   1.0, max =  5.0, N = 5, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :EFACs_EQUADs_max_objective, :F0, :F1, :F2, :DM, :PMRA, :PMDEC, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true, max=100.0),
    DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_3sigma], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

for i in 0:3
    if i == 0
        run_tempo_general(tf)
    else
        run_tempo_general(tf, just_refine=true)
    end

    chisqr_min = tf.grid.min[:chisqr]
    rc("mathtext",fontset="cm")
    rc("font", family="serif", size=12)
    fig, ax = subplots()
    # pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), rasterized=true)
    # cbar = colorbar(pclm)
    y_delta = tf.grid.y.values[2] - tf.grid.y.values[1]
    x_delta = tf.grid.x.values[2] - tf.grid.x.values[1]
    imsh = ax.imshow((tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), extent=(tf.grid.y.min - 0.5*y_delta, tf.grid.y.max + 0.5*y_delta, tf.grid.x.min - 0.5*x_delta, tf.grid.x.max + 0.5*x_delta), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), origin="lower", aspect="auto")
    cbar = colorbar(imsh)
    cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
    plot([], [], label=L"\Delta\chi^{2} (1\sigma, 2\sigma, 3\sigma)", "-",  color="red")
    cs2 = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:EFACs_EQUADs_max_objective], levels=[0.5, 0.75, 1.0, 1.25], linestyles=["-", "--", "-.", ":"], colors="green")
    ax.clabel(cs2, inline=true, fontsize=12, fmt="%.2f")
    plot([], [], label="AD statistics", "-",  color="green")
    annot_text = L"TNRedC = $%$TNRedC,\; \chi^2_\mathrm{min} = %$chisqr_min$"
    # Добавление текста с LaTeX
    text(
        0.05, 0.05, annot_text, transform=gca().transAxes,
        fontsize=12,
        bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    )
    ax.set_ylabel(tf.test_params.x.name, size=14)
    ax.set_xlabel(tf.test_params.y.name, size=14)
    title("Red noise + white noise")
    cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$", fontsize=14)
    legend(fontsize=12)
    tight_layout()
    savefig("saves/GR_PM_PL_RN$(TNRedC)_200d_$i.pdf", format="pdf")

    save("saves/tf_GR_PM_PL_RN$(TNRedC)_200d_$i.jld", "tf", tf)

end

# save("saves/tf_GR_RN_1y.jld", "tf", tf)
# tf = load("saves/tf_GR_RN30_200d_3.jld", "tf")

# tf = load("saves/tf_GR_RN84_200d_1.jld", "tf")

#-------------------------------------------------------------------------------------