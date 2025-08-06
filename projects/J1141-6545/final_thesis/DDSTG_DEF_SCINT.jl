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
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/ALPHA0_BETA0/MPA1",
    version = Tempo2(),
    par_file_init = "DDSTG_DEF_MPA1_SCINT_RN84.par",
    tim_file = "J1141-6545_pn_new.tim",
    flags = "-newpar -writeres -residuals",
    tparams = [TP("NITS", 4), TP("EOS", "MPA1")],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
)

general_settings = GravityToolsLight.GeneralTempoSettings(basic_settings)

test_params = TestParameters(
    Var(name = "STG_ALPHA0", min = -1e-4, max = -1e-1, N = 5, range_rule=:log),
    Var(name = "STG_BETA0",  min =  -6.0, max =  6.0,  N = 5, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

ref_sets = RefinementSettings(
    params_to_save = (:chisqr, :rms_post_fit, :rms_tn_post_fit, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT, :STG_ALPHA0, :STG_BETA0, :STG_ALPHA_P, :STG_BETA_P, :STG_K_P, :PBDOT_PHI_DIP, :PBDOT_G_QUAD),
    desired_refinement_level = 0,
    parallel = true,
    # LocalMinimaUnit(:chisqr_full, from_min=true, max=300.0),
    DiffContourUnit(:chisqr, diffs = [10.0], contours = [lvl_2sigma], from_min=true)
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
    x_log10_values = log10.(-tf.grid.x.values)
    fig, ax = subplots()
    y_delta = tf.grid.y.values[2] - tf.grid.y.values[1]
    x_delta = x_log10_values[2] - x_log10_values[1]
    imsh = ax.imshow((tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr]), extent=(tf.grid.y.min - 0.5*y_delta, tf.grid.y.max + 0.5*y_delta, x_log10_values[1] - 0.5*x_delta, x_log10_values[end] + 0.5*x_delta), cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20.0), origin="lower", aspect="auto")
    cbar = colorbar(imsh)
    cs = ax.contour(tf.grid.y.values, x_log10_values, tf.grid.vars[:chisqr] .- tf.grid.min[:chisqr], levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
    plot([], [], label=L"\Delta\chi^{2} (1\sigma, 2\sigma, 3\sigma)", "-",  color="red")
    # annot_text = L"\chi^2_\mathrm{min} = %$chisqr_min$"
    # Добавление текста с LaTeX
    # text(
    #     0.05, 0.05, annot_text, transform=gca().transAxes,
    #     fontsize=12,
    #     bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    # )
    ax.set_ylabel(L"\log_{10}|\alpha_0|", size=14)
    ax.set_xlabel(L"\beta_0", size=14)
    # title("Red noise + white noise")
    cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$", fontsize=14)
    legend(fontsize=12)
    tight_layout()
    savefig("saves/DDSTG_DEF_MPA1_SCINT_RN84_$i.pdf", format="pdf")

    save("saves/tf_DDSTG_DEF_MPA1_SCINT_RN84_$i.jld", "tf", tf)


end

# save("saves/tf_GR_RN_1y.jld", "tf", tf)
# tf = load("saves/tf_GR_RN30_200d_3.jld", "tf")

# tf = load("saves/tf_DDSTG_GR_SCINT_RN$(TNRedC)_200d_1.jld", "tf")

#-------------------------------------------------------------------------------------