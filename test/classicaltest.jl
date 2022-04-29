using GravityTools
using Measurements
using PyPlot
pygui(true)

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

K_params_obs = (Pb = 0.1022515592973, T0=55700.0, e0 = 0.087777023, omega0 = 204.753686, x0 = 1.415028603)                

PK_params_obs = (k = (16.899323 ± 0.000013) / 360 * K_params_obs.Pb/365.25,
                gamma = 0.384045e-3 ± 0.000094e-3,
                Pbdot = -1.247752e-12 ± 0.000079e-12,
                r = 6.162e-6 ± 0.021e-6,
                s = 0.999936 ± 0.000010)

ct = ClassicalTest("Double Pulsar", K_params_obs, PK_params_obs)

X_params_obs = (m2 = 1.25 ± 0.01, q = 1.07 ± 1.01)
ct = ClassicalTest("Double Pulsar", K_params_obs, PK_params_obs, X_params_obs)

#-------------------------------------------------------------------------------------

theory = DEF()
eosname = :MPA1
bnsys = BinarySystem(:NS, :NS)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)
pf.theory.alpha0 = -0.0000
pf.theory.beta0  = -0.0000
interpolate_mgrid!(pf)

ct = ct_dataset["Double Pulsar"]
ct = ct_dataset["J2222−0137"]
pf.bnsys.K_params = ct.K_params_obs

pf.bnsys.psr.mass = 1.338185
pf.bnsys.comp.mass = 1.248868
interpolate_bnsys!(pf)
GravityTools.get_chisqr(ct, pf)
check_terms_in_chisqr(ct, pf)

#pf.bnsys.psr.mass = 1.4
#pf.bnsys.comp.mass = 1.3
#interpolate_bnsys!(pf)
#GravityTools.get_chisqr(ct, pf)

find_initial_masses(ct, pf)
GravityTools.get_chisqr(ct, pf)
check_terms_in_chisqr(ct, pf)

find_best_masses(ct, pf)
GravityTools.get_chisqr(ct, pf)
check_terms_in_chisqr(ct, pf)

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
ct.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
ct.N_refinement = 6
ct.CL = 0.90
calculate!(ct, pf)

#-------------------------------------------------------------------------------------
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)
ct = ct_dataset["J2222−0137"]
pf.bnsys.K_params = ct.K_params_obs

pf.bnsys.psr.mass = 1.820
pf.bnsys.comp.mass = 1.3153
interpolate_bnsys!(pf)
GravityTools.get_chisqr(ct, pf)
check_terms_in_chisqr(ct, pf)

find_best_masses(ct, pf)
GravityTools.get_chisqr(ct, pf)
check_terms_in_chisqr(ct, pf)

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
ct.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
ct.N_refinement = 6
ct.CL = 0.90
calculate!(ct, pf)

#-------------------------------------------------------------------------------------

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(ct.grid.y, ct.grid.x, ct.grid.value[:chisqr] .- ct.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20), rasterized=true)
cs = ax.contour(ct.grid.y, ct.grid.x, ct.grid.value[:chisqr] .- ct.grid.params[:chisqr_min], levels=[4.60], colors="black")
#clabel(cs, cs.levels, fmt=Dict(lvl_2σ => "2σ"))
#cs.collections[1].set_label("$(ct.name)")
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
axhline(y=log10(3.4e-3), color="gray", label="Cassini")
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}(|\alpha_0|)", size=16)
legend(fontsize=12)
ylim(-4.0,-1.0)
tight_layout()

plot([], [], label="PK test", color="black")
grid_ddstg = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")
cs2 = ax.contour(grid_ddstg.y, grid_ddstg.x, grid_ddstg.value[:chisqr] .- grid_ddstg.params[:chisqr_min], levels=[4.60], colors="red")
plot([], [], label="DDSTG test", color="red")
legend(fontsize=12)

#fig, ax = subplots()
#pclm = ax.scatter(ct.grid.value[:m1] .+ ct.grid.value[:m2] .- 2.587053, ct.grid.value[:m1] .- ct.grid.value[:m2] .- 0.089317, alpha=0.4, c=ct.grid.value[:chisqr].-ct.grid.params[:chisqr_min], norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20), linewidth=0.7)
#cbar = colorbar(pclm)
#cbar.set_label(L"$\Delta\chi^{2}$")
#xlabel(L"m_{\mathrm{p}},\, [M_{\odot}]", size=16)
#ylabel(L"m_{\mathrm{c}},\, [M_{\odot}]", size=16)
#tight_layout()
#
#
#fig, ax = subplots()
#pclm = ax.scatter(ct.grid.value[:m1] .+ ct.grid.value[:m2], ct.grid.value[:chisqr], alpha=0.4, c=ct.grid.value[:chisqr].-ct.grid.params[:chisqr_min], norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20), linewidth=0.7)
#cbar = colorbar(pclm)
#cbar.set_label(L"$\Delta\chi^{2}$")
#xlabel(L"m_{\mathrm{p}},\, [M_{\odot}]", size=16)
#ylabel(L"m_{\mathrm{c}},\, [M_{\odot}]", size=16)
#tight_layout()


