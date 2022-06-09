using Revise
using GravityTools
using Measurements
using PyPlot
using JLD
pygui(true)

eosname = :WFF1
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

pf.eosname = :WFF1
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