using GravityTools
using PyPlot
pygui(true)

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)
K_params_obs = ct_dataset["J2222−0137"].K_params_obs
PK_params_obs = ct_dataset["J2222−0137"].PK_params_obs
pf.bnsys.K_params = K_params_obs

pf.bnsys.psr.mass = 1.820
pf.bnsys.comp.mass = 1.3153
interpolate_bnsys!(pf)

m1_arr = collect(LinRange(1.0, 1.5, 11))
mm = MMDiagram(m1_arr, K_params_obs, PK_params_obs)

calculate!(mm, pf)