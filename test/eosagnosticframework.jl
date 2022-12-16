using Revise
using Distributed
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics

pygui(true)

addprocs(4)

@everywhere using GravityTools 


#-------------------------------------------------------------------------------------

cd("/Users/abatrakov/Documents/PhD_work/projects/test_tempo2/")

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 2),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 2)
    )

tsets = TempoSettings(
    par_file_init = "test_fake.par",
    tim_file = "test_fake.simulate",
    add_flag = "-newpar",
    fit_XPBDOT = false,
    nits_first_step = 3,
    gain_fisrt_step = 1.0
    )

gsets = GridSetttings(
    N_refinement = 6,
    CL = [0.90],
    refinement_type = "contour",
    delta_chisqr_max = 10,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

#theory = DEF(0.0, 0.0)