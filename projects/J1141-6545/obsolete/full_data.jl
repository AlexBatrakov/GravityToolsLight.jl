using Contour
using JLD
using ColorSchemes
using Statistics
using DelimitedFiles
using Distributions
using StructArrays
using StatsBase
using Interpolations
using QuadGK
using Optim
using Roots
using HypothesisTests
using KernelDensity
using QuadGK

using Revise
using GravityToolsLight
using PyPlot
pygui(true)

#------------------------------------------------------------------------------------- 
# WE HAVE PRIOR INFORMATION ABOUT XPBDOT
# use sweep or global iterations over XPBDOT to check likelihood
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# FIRSTLY WE FIT FOR XPBDOT IN TEMPO
#-------------------------------------------------------------------------------------
basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_data",
    version = Tempo2(),
    par_file_init = "DDSTG_best.par",
    tim_file = "J1141-6545_pn.tim",
    flags = "-nobs 23000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 1), TP("XPBDOT", flag=1), TP("M2", flag=1)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )
    
results_basic = run_tempo_basic(basic_settings)

# BEST FIT 
chisqr_min = results_basic.last_internal_iteration.result.chisqr
xpbdot_mu = results_basic.last_internal_iteration.result.XPBDOT.post_fit
xpbdot_sigma = results_basic.last_internal_iteration.result.XPBDOT.uncertainty

# WE OBTAINED MU AND SIGMA FOR NORMAL LIKELIHOOD
# chisqr_min, xpbdot_mu, xpbdot_sigma = 20914.19, -8.29481968690372e-15,  4.7427e-15
# BUILD NORMAL LIKELIHOOD
xpbdot_normal_likelihood(xpbdot) = pdf(Normal(xpbdot_mu, xpbdot_sigma), xpbdot)
# BUILD NORMAL LIKELIHOOD IN CHISQR REPRESENTATION
delta_chisqr_normal_likelihood(xpbdot) = -2 * log(xpbdot_normal_likelihood(xpbdot)) + 2 * log(xpbdot_normal_likelihood(xpbdot_mu))