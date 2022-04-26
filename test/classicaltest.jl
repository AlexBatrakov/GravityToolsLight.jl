using GravityTools

eosname = :MPA1
path_to_grids = "/Users/abatrakov/Documents/PhD_work/projects/computed_grids"

grid = read_DEFGrid(eosname, path_to_grids)
mgrid = interpolate_DEFMassGrid(grid, -0.00123, -4.123)
interpolate_NS(mgrid, 1.5)

#-------------------------------------------------------------------------------------

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


K_params = (Pb = 0.10225156248, T0=58002.0192787585, e0 = 0.0877775, omega0 = 88.69, x0 = 1.415032)
PK_params_obs   = (k = 16.89947 / 360 * K_params.Pb/365.25, gamma = 0.3856e-3, Pbdot = -1.252e-12, r = 6.21e-6, s = 0.99974)
PK_params_delta = (k =  0.00068 / 360 * K_params.Pb/365.25, gamma = 0.0026e-3, Pbdot =  0.017e-12, r = 0.33e-6, s = 0.00039)

cl = ClassicalTest(K_params, PK_params_obs, PK_params_delta)

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
cl.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
cl.N_refinement = 6
cl.CL = 0.90

pf.theory.alpha0 = -0.0001
pf.theory.beta0  = -0.0000
pf.bnsys.psr.mass = 1.338185
pf.bnsys.comp.mass = 1.248868
interpolate_bnsys!(pf)

pf.bnsys.psr.mass = 1.4
pf.bnsys.comp.mass = 1.3

find_initial_masses(cl, pf)
