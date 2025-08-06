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


using Distributed
addprocs(8)

@everywhere using Revise
@everywhere using GravityToolsLight

#----------------------------------------------------------------------------------

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data_thesis/DDSTG_PHEN/PREC_ABE",
    version = Tempo2(),
    par_file_init = "DDSTG_PHEN_PREC_ABE.par",
    tim_file = "J1141-6545_pn_clean.tim",
    flags = "-newpar",
    tparams = [TP("GAIN", 0.2), TP("NITS", 10), TP("GAIN_ONE_AFTER_NITS", 5)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "J_DELTA",     min = 0.0,  max = 0.6,  N = 7, range_rule=:lin),
    Var(name = "PREC_PHIDOT", min = -5.0, max = 15.0, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :OM2DOT, :M2, :GAMMA, :I, :IDOT, :I2DOT, :J_I, :J_DELTA, :PREC_OM_SEC, :PREC_OMDOT_SEC, :PREC_PHI, :PREC_PHIDOT, :P_LAMBDA, :P_ETA, :P_DELTA, :P_PHI, :P_PHIDOT, :ABE_X, :ABE_XDOT, :ABE_X2DOT, :ABE_PBDOT, :ABE_E, :ABE_EDOT),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true, max=20.0)
    # DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_3sigma], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

for i in 0:6
    if i == 0
        run_tempo_general(tf)
    else
        run_tempo_general(tf, just_refine=true)
    end

    try
    chisqr_min = tf.grid.min[:chisqr]
    rc("mathtext",fontset="cm")
    rc("font", family="serif", size=12)
    fig, ax = subplots()
    # pclm = ax.pcolormesh(tf.grid.y.values, tf.grid.x.values, (tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=100), rasterized=true)
    # cbar = colorbar(pclm)
    y_delta = tf.grid.y.values[2] - tf.grid.y.values[1]
    x_delta = tf.grid.x.values[2] - tf.grid.x.values[1]
    imsh = ax.imshow((tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), extent=(tf.grid.y.min - 0.5*y_delta, tf.grid.y.max + 0.5*y_delta, tf.grid.x.min - 0.5*x_delta, tf.grid.x.max + 0.5*x_delta), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20), origin="lower", aspect="auto")
    cbar = colorbar(imsh)
    cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
    plot([], [], label=L"\Delta\chi^{2}", "-",  color="red")
    cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:M2], levels=[0.95, 1.0, 1.05], linestyles=["-.", "-", "--"], colors="green")
    plot([], [], label=L"m_c", "-",  color="green")
    cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:GAMMA], levels=[0.0007, 0.00075, 0.0008], linestyles=["-.", "-", "--"], colors="black")
    plot([], [], label=L"\gamma", "-",  color="black")
    annot_text = L"$\chi^2_\mathrm{min} = %$chisqr_min$"
    # Добавление текста с LaTeX
    text(
        0.05, 0.05, annot_text, transform=gca().transAxes,
        fontsize=12,
        bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    )
    ax.set_ylabel(tf.test_params.x.name, size=14)
    ax.set_xlabel(tf.test_params.y.name, size=14)
    title("Phenomenological precession and aberration")
    cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$", fontsize=14)
    legend(fontsize=12)
    tight_layout()
    savefig("saves/DDSTG_PHEN_PREC_ABE_$i.pdf", format="pdf")

    catch e

    end
    save("saves/DDSTG_PHEN_PREC_ABE_$i.jld", "tf", tf)
    
end

# tf = load("saves/DDSTG_PHEN_PREC_ABE_$i.jld", "tf")