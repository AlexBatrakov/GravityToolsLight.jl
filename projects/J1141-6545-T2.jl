using Revise
using GravityTools
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 2),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 2)
    )

tsets = TempoSettings(
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    add_flag = "-nobs 22000 -newpar",
    fit_XPBDOT = false,
    nits_first_step = 3,
    gain_fisrt_step = 1.0
    )

gsets = GridSetttings(
    N_refinement = 3,
    CL = [0.90],
    refinement_type = "contour",
    delta_chisqr_max = 10,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Science/Software/tempo2/T2runtime/ddstg_data")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

#obs_params = obs_params_dataset["J2222−0137_sim_DD"]
#pf.bnsys.K_params = obs_params.K
#pf.bnsys.psr.mass = obs_params.masses_init.m1
#pf.bnsys.comp.mass = obs_params.masses_init.m2
#interpolate_bnsys!(pf)

calculate_t2!(tf)