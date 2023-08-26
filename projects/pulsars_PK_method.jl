#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :BSk22
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Work/PhD/computed_grids_fine")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1738+0333_Guo"]
#obs_params = obs_params_dataset["J2222−0137_Guo_DDK"]
#obs_params = obs_params_dataset["Triple System"]
#obs_params = obs_params_dataset["Double Pulsar"]
#obs_params = obs_params_dataset["J1141-6545_DDFWHE"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = obs_params.masses_init.m1
pf.bnsys.comp.mass = obs_params.masses_init.m2
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

##find_initial_masses(obs_params, pf)
#GravityTools.get_chisqr(obs_params, pf)
#check_terms_in_chisqr(obs_params, pf)

pf.theory.alpha0 = -0.0
pf.theory.beta0  = -0.0
interpolate_mgrid!(pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

optimize_PK_method(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "BSk22",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
    )

gsets = GridSetttings(
    N_refinement = 6,
    CL = [0.90],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_J1738+0333_Guo_BSk22.jld",  "grid", pkf.grid)
#save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Double_Pulsar_BSk22.jld",  "grid", pkf.grid)
#save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Triple_System_BSk22.jld",  "grid", pkf.grid)
#save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDFWHE_J1141-6545_BSk22.jld",  "grid", pkf.grid)