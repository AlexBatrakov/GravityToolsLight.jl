using GravityTools

eosname = :MPA1
path_to_grids = "/Users/abatrakov/Documents/PhD_work/projects/computed_grids"

grid = read_DEFGrid(eosname, path_to_grids)
mgrid = interpolate_DEFMassGrid(grid, -0.00123, -4.123)
interpolate_NS(mgrid, 1.5)

#-------------------------------------------------------------------------------------
# check PhysicalFramework

theory = DEF()
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

pf.eosname = :MPA1
read_grid!(pf)

pf.theory.alpha0 = -0.0123
pf.theory.beta0  = -4.1234

interpolate_mgrid!(pf)

pf.bnsys.psr.mass = 1.5
pf.bnsys.comp.mass = 1.3

interpolate_psr!(pf)
interpolate_comp!(pf)

pf.bnsys.name = "Double Pulsar"
pf.bnsys.K_params = (Pb = 0.10225156248, T0=58002.0192787585, e0 = 0.0877775, omega0 = 88.69, x0 = 1.415032)

calculate_PK_params!(pf)

pf.bnsys.comp.mass = 1.3
interpolate_bnsys!(pf)

#-------------------------------------------------------------------------------------
theory = DEF()
eosname = :MPA1
bnsys = BinarySystem(:NS, :NS)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)
pf.theory.alpha0 = -0.0001
pf.theory.beta0  = -4.4000
interpolate_mgrid!(pf)


K_params_obs = (Pb = 0.10225156248, T0=58002.0192787585, e0 = 0.0877775, omega0 = 88.69, x0 = 1.415032)
#PK_params_exact   = (k = 16.89947 / 360 * K_params.Pb/365.25, gamma = 0.3856e-3, Pbdot = -1.252e-12, r = 6.21e-6, s = 0.99974)
#PK_params_delta = (k =  0.00068 / 360 * K_params.Pb/365.25, gamma = 0.0026e-3, Pbdot =  0.017e-12, r = 0.33e-6, s = 0.00039)

PK_params_obs = (k = (16.89947 ± 0.00068) / 360 * K_params_obs.Pb/365.25,
                gamma = 0.3856e-3 ± 0.0026e-3,
                Pbdot = -1.252e-12 ± 0.017e-12,
                r = 6.21e-6 ± 0.33e-6,
                s = 0.99974 ± 0.00039)

K_params_obs = (Pb = 0.1022515592973, T0=55700.0, e0 = 0.087777023, omega0 = 204.753686, x0 = 1.415028603)                

PK_params_obs = (k = (16.899323 ± 0.000013) / 360 * K_params_obs.Pb/365.25,
                gamma = 0.384045e-3 ± 0.000094e-3,
                #Pbdot = -1.247920e-12 ± 0.000078e-12,
                Pbdot = -1.247752e-12 ± 0.000079e-12,
                r = 6.162e-6 ± 0.021e-6,
                s = 0.999936 ± 0.000010)

cl = ClassicalTest(K_params_obs, PK_params_obs)

#X_params_obs = (m2 = 1.25 ± 0.01, q = 1.07 ± 1.01)
#cl = ClassicalTest(K_params_obs, PK_params_obs, X_params_obs)

pf.bnsys.name = "Double Pulsar"
pf.bnsys.K_params = K_params_obs

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
cl.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
cl.N_refinement = 6
cl.CL = 0.90

pf.bnsys.psr.mass = 1.338185
pf.bnsys.comp.mass = 1.248868
interpolate_bnsys!(pf)
GravityTools.get_chisqr(cl, pf)
check_terms_in_chisqr(cl, pf)

#pf.bnsys.psr.mass = 1.4
#pf.bnsys.comp.mass = 1.3
#interpolate_bnsys!(pf)
#GravityTools.get_chisqr(cl, pf)

find_initial_masses(cl, pf)
GravityTools.get_chisqr(cl, pf)

find_best_masses(cl, pf)
GravityTools.get_chisqr(cl, pf)
check_terms_in_chisqr(cl, pf)