using Revise
using Distributed
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics

pygui(true)

addprocs(8)

@everywhere using GravityTools 


#-------------------------------------------------------------------------------------

cd("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/science_paper")

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -2.0, N = 2),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 2)
    )

tsets = TempoSettings(
    par_file_init = "J1141-6545_until_2018.par",
    tim_file = "J1141-6545_until_2018.tim",
    add_flag = "-newpar -nobs 22000",
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

#theory = DEF(0.0, 0.0)#eosname = :MPA1
#bnsys = BinarySystem(:NS, :WD)
#sets = Settings(ENV["TEMPO2"] * "/data_ddstg")
#pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
#read_grid!(pf)

eos_list = ["SLy", "APR3", "APR4", "WFF1", "WFF2", "ENG", "MPA1", "MS1", "MS1b", "H4", "ALF2"]
eos_agn_test = EOSAgnosticTest(eos_list, tf)
calculate!(eos_agn_test)
#plot_test(eos_agn_test); savefig("refinement=$(eos_agn_test.tf_array[1].grid.N_refinement)")
#calculate!(eos_agn_test, add_refinement=2)