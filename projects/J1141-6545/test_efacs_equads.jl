using Revise
using GravityToolsLight
using PyPlot
pygui(true)


basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

results_basic = run_tempo_basic(basic_settings)

EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings)


res = residuals[indices]
unc = uncertainties[indices]

res_std = res ./ std(res)
unc_std = unc ./ std(unc)

trans(unc, EFAC, EQUAD) = transform_uncertainty.(unc, EFAC, EQUAD) ./ std(transform_uncertainty.(unc, EFAC, EQUAD))


plot(res_std, res_std ./ unc_std, ".")
plot(res_std, res_std ./ trans(unc, 0.0, 100.0), ".")
plot(res_std, res_std ./ trans(unc, 4.0, 0.0), ".")


alpha(res_std, unc, EFAC, EQUAD) = atan.(res_std ./ trans(unc, EFAC, EQUAD) ./     std(residuals[indices] ./ uncertainties[indices]) ./ uncertainties[indices])

plot(atan.(std(residuals[indices] ./ uncertainties[indices]) ./ uncertainties[indices]), residuals[indices] ./ uncertainties[indices] ./ std(residuals[indices] ./ uncertainties[indices]), ".")


plot(unc_std, res_std, ".")

plot(res_std ./ unc_std, atan.(res_std ./ unc_std, res_std),".")