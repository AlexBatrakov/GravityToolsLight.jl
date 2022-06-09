using Revise
using GravityTools
using PyPlot
using JLD
pygui(true)

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1141-6545_DD"]
pf.bnsys.K_params = obs_params.K

pf.bnsys.psr.mass = 1.273322
pf.bnsys.comp.mass = 1.016393
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_initial_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
    )

gsets = GridSetttings(
    N_refinement = 1,
    CL = 0.90,
    plot_type = "nice"
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE"]
pf.bnsys.K_params = obs_params.K

pf.bnsys.psr.mass = 1.241308
pf.bnsys.comp.mass = 0.950913
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

pf.bnsys.psr.mass =  0.9819032
pf.bnsys.comp.mass =0.8311580
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

##find_initial_masses(obs_params, pf)
#GravityTools.get_chisqr(obs_params, pf)
#check_terms_in_chisqr(obs_params, pf)

pf.theory.alpha0 = -0.0001
pf.theory.beta0  = -6.0
interpolate_mgrid!(pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

optimize_PK_method(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
    )

gsets = GridSetttings(
    N_refinement = 6,
    CL = 0.90,
    plot_type = "nice"
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)
#-------------------------------------------------------------------------------------

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1141-6545_DD"]
pf.bnsys.K_params = obs_params.K

pf.bnsys.psr.mass = 1.273322
pf.bnsys.comp.mass = 1.016393
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_initial_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "MPA1",
    alpha0 = 0.0,
    beta0 = 0.0,
    param2 = (name = "Pbdot", min = -0.32e-12, max = -0.45e-12, N = 9),
    param1 = (name = "gamma", min = 0.0004, max = 0.0011, N = 9)
    )

gsets = GridSetttings(
    N_refinement = 1,
    CL = 0.90,
    plot_type = "nice"
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)
#-------------------------------------------------------------------------------------

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_5σ), rasterized=true)
cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="black")
clabel(cs2, cs2.levels)
plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

cs = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="red")
clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")
ax.set_ylabel(get_label(pkf.test.param1.name), size=16)
ax.set_xlabel(get_label(pkf.test.param2.name), size=16)
ax.axhline(y=log10(3.4e-3), color="gray", label="Cassini", ls="--", zorder=5)

title(L"\mathrm{PK\ method}")
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()