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

basic_settings = settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR/PREC",
    version = Tempo2(),
    par_file_init = "DDSTG_GR_SCINT_DM_PREC.par",
    tim_file = "J1141-6545_pn_new.tim",
    flags = "-newpar",
    tparams = [TP("NITS", 10)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 3,
    nits = [5,10,1],
    gain = [0.2, 0.5, 1.0],
    print_output = true
)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "PREC_XI", min = 0.0, max = 360.0, N = 7, range_rule=:lin),
    Var(name = "C_CHI",   min = 0.0, max = 3.0, N = 7, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :pre_post, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :DM, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :X2DOT, :OMDOT, :OM2DOT, :M2, :MTOT, :GAMMA, :I, :IDOT, :I2DOT, :J_I, :J_DELTA, :OMDOT_SQ, :OMDOT_PN, :PREC_PHI, :PREC_XI, :PREC_PHIDOT, :C_CHI, :C_DELTA),
    desired_refinement_level = 0,
    parallel = true,
    # FullUnit(:chisqr)
    LocalMinimaUnit(:chisqr, from_min=true, max=20.0),
    DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_1sigma], from_min=true)
    )

tf = GeneralTempoFramework(general_settings, test_params, ref_sets)

for i in 0:6
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
    imsh = ax.imshow((tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), extent=(tf.grid.y.min - 0.5*y_delta, tf.grid.y.max + 0.5*y_delta, tf.grid.x.min - 0.5*x_delta, tf.grid.x.max + 0.5*x_delta), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20.0), origin="lower", aspect="auto")
    cbar = colorbar(imsh)
    cs = ax.contour(tf.grid.y.values, tf.grid.x.values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
    plot([], [], label=L"\Delta\chi^{2} (1\sigma, 2\sigma, 3\sigma)", "-",  color="red")
    annot_text = L"$\chi^2_\mathrm{min} = %$chisqr_min$"
    # Добавление текста с LaTeX
    text(
        0.05, 0.05, annot_text, transform=gca().transAxes,
        fontsize=12,
        bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    )
    ax.set_ylabel(tf.test_params.x.name, size=14)
    ax.set_xlabel(tf.test_params.y.name, size=14)
    title("Precession")
    cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$", fontsize=14)
    legend(fontsize=12)
    tight_layout()
    savefig("saves/DDSTG_GR_SCINT_DM_PREC_$i.pdf", format="pdf")
    save("saves/DDSTG_GR_SCINT_DM_PREC_$i.jld", "tf", tf)
    
end