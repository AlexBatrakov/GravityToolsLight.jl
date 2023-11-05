using Revise
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
using DelimitedFiles
using Distributions
#pygui(true)

using Revise
using GravityToolsLight
using PyPlot
pygui(true)

# using Distributed
# addprocs(8)
# @everywhere using GravityToolsLight
# #pkill -f /opt/julia-1.8.2/bin/julia

#-------------------------------------------------------------------------------------

# tsets = GeneralTempoSettings(
#     work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/model_WN",
#     version = Tempo2(),
#     par_file_init = "J1141-6545_full_DDSTG_RN.par",
#     tim_file = "J1141-6545_full.tim",
#     flags = "-nobs 34000 -newpar -writeres -residuals",
#     tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
#     keys = TempoKeys(silent=true, print_output=true, iterative_mode=true, fit_EFACs_EQUADs=true),
#     iters = 3,
#     nits = [3,3,3],
#     tparams_local = [
#         [TP("M2", 1.0, flag=0)],
#         [TP("M2", flag=1)],
#         [TP("M2", flag=1)]
#         ]
#     )

# parsed_results, output, stderr_output = run_tempo_single(tsets)

# GravityToolsLight.plot_iterations_data(parsed_results)

# plot_fit_results(parsed_results)

basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/model_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_full_DDSTG_RN.par",
    tim_file = "J1141-6545_full.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = TempoKeys(silent=true, print_output=true, iterative_mode=true, fit_EFACs_EQUADs=true)
    )

parsed_output, output, stderr_output = run_tempo_basic(basic_settings)
plot_fit_results(parsed_results)

global_iter_settings = GlobalIterationSettings(
    iters = 3,
    nits = [3,3,3],
    tparams_local = [
        [TP("MTOT", flag=0), TP("M2", 1.0, flag=0)],
        [TP("MTOT", flag=1), TP("M2", 1.0, flag=0)],
        [TP("M2", flag=1)]
        ]
    )

results = run_tempo_global_iter(basic_settings, global_iter_settings)
plot_fit_results(results)

parameter_sweep_settings = ParameterSweepSettings(
    "XPBDOT",
    [-1e-14, 0.0, 1e-14]
)

general_settings = NewGeneralTempoSettings(basic_settings, global_iter_settings, parameter_sweep_settings)

# test = GeneralTest(
#     psrname = "J1141-6545",
#     eosname = "BSk22",
#     param1 = (name = "TNRedAmp", min = 8.0, max = 12.0, N = 5),
#     param2 = (name = "TNRedGam", min = 1.0, max = 4.0, N = 4)
#     )

test_params = TestParameters(
    Var(name = "TNRedAmp", min = 8.0, max = 12.0, N = 5, range_rule=:lin),
    Var(name = "TNRedGam", min = 1.0, max = 4.0, N = 4, range_rule=:lin),
    ValueVariable[],
    RangeVariable[]
)

# gsets = GridSetttings(
#     N_refinement = 0,
#     CL = [0.90],
#     refinement_type = "full",
#     delta_chisqr_max = 10,
#     delta_chisqr_diff = 1.0,
#     gr_in_chisqr = true
#     )

ref_sets = RefinementSettings(
    desired_refinement_level = 2,
    parallel = false,
    FullUnit(:chisqr)
#    DiffUnit(:val1, diff = 0.1),
#    ContourUnit(:val1, contours = [0.5])
#    DiffContourUnit(:val1, diff = 0.05, contours = [0.5])
    )

tf = GeneralTempoFramework(tsets, test_params, ref_sets)

calculate!(tf)



#-------------------------------------------------------------------------------------

par_file = TempoParFile(tf.tsets.par_file_init)
par_file.tparams[:TNRedAmp].value = -8.0
par_file.tparams[:TNRedGam].value = 3.0
tparam = GeneralTempoParameter("TNEF -be CASPSR", 0.991577)
extend_par_file!(par_file, tparam)
write_par_file(par_file, "test.par")

#-------------------------------------------------------------------------------------