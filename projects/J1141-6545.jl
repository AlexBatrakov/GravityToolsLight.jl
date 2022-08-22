using Revise
using GravityTools
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J1141-6545",
    alpha0 = -0.0,
    beta0 =  -0.0,
#    param1 = (name = "M2", min = 0.0, max = 2.0, N = 9),
    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
#    param2 = (name = "DTHETA", min = -500.0, max = 500.0, N = 9),
#    param1 = (name = "OMDOT", min = 5.3098, max = 5.3115, N = 9),
#    param2 = (name = "COSI", min = 0.0, max = 1.0, N = 9),
#    param1 = (name = "H3", min = 0.0, max = 0.000006, N = 9),
#    param1 = (name = "VARSIGMA", min = 0.0, max = 1.5, N = 9),
#    param2 = (name = "XDOT", min = -1.0, max = 1.5, N = 9)
    )

tsets = TempoSettings(
#    par_file_init = "J1141-6545_T1_DD.par",
    par_file_init = "J1141-6545_T1_DDFWHE.par",
    tim_file = "J1141-6545_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = false,
    nits_first_step = 5,
    gain_fisrt_step = 1.0,
    )

gsets = GridSetttings(
    N_refinement = 4,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = false
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

save("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DD_GAMMA_PBDOT.jld", "grid", tf.grid)

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J1141-6545",
    alpha0 = -0.0025,
    beta0 =  6.0,
    param1 = (name = "M2", min = 1.016393 -3*0.024070, max = 1.016393 + 3*0.024070, N = 9),
    param2 = (name = "MTOT", min = 2.289715 - 3*0.000064, max = 2.289715 + 3*0.000064, N = 9),
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
#    param2 = (name = "DTHETA", min = -500.0, max = 500.0, N = 9),
#    param1 = (name = "OMDOT", min = 5.3098, max = 5.3115, N = 9),
#    param1 = (name = "COSI", min = 0.0, max = 1.0, N = 9),
#    param1 = (name = "H3", min = 0.0, max = 0.000006, N = 9),
#    param1 = (name = "VARSIGMA", min = 0.0, max = 1.5, N = 9),
#    param2 = (name = "XDOT", min = -1.0, max = 1.5, N = 9)
#    param1 = (name = "M2", min = 0.0, max = 2.0, N = 9),
#    param2 = (name = "COSI", min = 0.0, max = 1.0, N = 9)
    )

tsets = TempoSettings(
#    par_file_init = "J1141-6545_T1_DD.par",
    par_file_init = "J1141-6545_T1_DDSTG.par",
    tim_file = "J1141-6545_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = false,
    nits_first_step = 5,
    gain_fisrt_step = 1.0,
    )

gsets = GridSetttings(
    N_refinement = 4,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = false
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

save("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DDSTG_MM.jld", "grid", tf.grid)

#-------------------------------------------------------------------------------------

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1141-6545_DDFWHE"]
pf.bnsys.K_params = obs_params.K

pf.theory.alpha0 = -0.0
pf.theory.beta0  = 0.0
interpolate_mgrid!(pf)

pf.bnsys.psr.mass = 1.273322
pf.bnsys.comp.mass = 1.016393
interpolate_bnsys!(pf)
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
#    alpha0 = -0.0025,
#    beta0 =  6.0,
#    param2 = (name = "Pbdot", min = -0.32e-12, max = -0.45e-12, N = 9),
#    param1 = (name = "gamma", min = 0.0004, max = 0.0011, N = 9),
#    param1 = (name = "m2", min = 0.1, max = 2.0, N = 9),
#    param2 = (name = "s", min = 0.1, max = 1.0, N = 9)
    )

gsets = GridSetttings(
    N_refinement = 6,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

save("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_PK_GAMMA_PBDOT.jld", "grid", pkf.grid)

#-------------------------------------------------------------------------------------\

grid_DD_GAMMA_PBDOT = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DD_GAMMA_PBDOT.jld", "grid")
grid_PK_GAMMA_PBDOT = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_PK_GAMMA_PBDOT.jld", "grid")
grid_PK_GAMMA_PBDOT_DEF = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_PK_GAMMA_PBDOT_DEF.jld", "grid")
grid_DDSTG_MM = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DDSTG_MM.jld", "grid")
grid_circle = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_circle.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
#pclm = ax.pcolormesh(grid_DD_GAMMA_PBDOT.y *1e-12, grid_DD_GAMMA_PBDOT.x, grid_DD_GAMMA_PBDOT.value[:chisqr_cut] .- grid_DD_GAMMA_PBDOT.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)

cs_DD_GAMMA_PBDOT = ax.contour(grid_DD_GAMMA_PBDOT.y *1e-12, grid_DD_GAMMA_PBDOT.x, grid_DD_GAMMA_PBDOT.value[:chisqr_cut] .- grid_DD_GAMMA_PBDOT.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(7)], linestyles=["-"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")

cs_circle = ax.contour(grid_circle.y, grid_circle.x, grid_circle.value[:chisqr] .- grid_circle.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(7)], linestyles=["--"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color=cm(7))


cs_PK_GAMMA_PBDOT = ax.contour(grid_PK_GAMMA_PBDOT.y, grid_PK_GAMMA_PBDOT.x, grid_PK_GAMMA_PBDOT.value[:chisqr] .- grid_PK_GAMMA_PBDOT.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(1)])
plot([], [], label=L"\mathrm{PK,\ GR}", color=cm(1))

#inds = filter(i -> (grid_DDSTG_MM.value[:chisqr_cut][i] .- grid_DDSTG_MM.params[:chisqr_min] <= lvl_68CL), 1:length(grid_DDSTG_MM.value[:chisqr]))
#pclm2 = ax.plot(grid_DDSTG_MM.value[:PBDOT][inds] .* 1e-12, grid_DDSTG_MM.value[:GAMMA][inds], color=cm(1), linewidth=1.0)
#plot([], [], label=L"\mathrm{DDSTG,\ GR}", color=cm(1))
pclm2 = ax.plot([-3.875734e-13, -3.850819e-13], [0.000790051, 0.0007574919], color=cm(2), ls=":", label=L"\mathrm{DDSTG,\ GR}", linewidth=2.0)

cs_PK_GAMMA_PBDOT_DEF = ax.contour(grid_PK_GAMMA_PBDOT_DEF.y, grid_PK_GAMMA_PBDOT_DEF.x, grid_PK_GAMMA_PBDOT_DEF.value[:chisqr] .- grid_PK_GAMMA_PBDOT_DEF.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(0)])
plot([], [], label=L"\mathrm{PK,\ DEF(-0.0025, 6.0)}", color=cm(0))

xlim(-4.3e-13, -3.4e-13)
ylim(0.0006, 0.00095)

ax.set_ylabel(get_label("GAMMA"), size=16)
ax.set_xlabel(get_label("Pbdot"), size=16)
legend(fontsize=12)
#ax.invert_xaxis()
tight_layout()


fig, ax = subplots()
pclm = ax.pcolormesh(grid_DDSTG_MM.y, grid_DDSTG_MM.x, grid_DDSTG_MM.value[:chisqr_cut] .- grid_DDSTG_MM.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)

fig, ax = subplots()
inds = filter(i -> (grid_DDSTG_MM.value[:chisqr_cut][i] .- grid_DDSTG_MM.params[:chisqr_min] <= lvl_68CL), 1:length(grid_DDSTG_MM.value[:chisqr]))
pclm = ax.plot(grid_DDSTG_MM.value[:PBDOT][inds], grid_DDSTG_MM.value[:GAMMA][inds], color="orange", ".")
#pclm = ax.scatter(grid_DDSTG_MM.value[:PBDOT][inds], grid_DDSTG_MM.value[:GAMMA][inds], alpha=0.1, linewidth=0.05, c=(grid_DDSTG_MM.value[:chisqr_cut][inds] .- grid_DDSTG_MM.params[:chisqr_min]), norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)
#cbar = colorbar(pclm)
#cbar.set_alpha(1)
#cbar.draw_all()
ax.set_yscale("symlog", linthresh=1e-5)
xlim(0.0,10.0)
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta \chi^2_{\mathrm{DDSTG}}", size=16)
ax.set_ylabel(L"\alpha_{A,\mathrm{PK}} - \alpha_{A,\mathrm{DDSTG}}", size=16)
legend(fontsize=12)
tight_layout()

#-------------------------------------------------------------------------------------\

grid_DD_M2_COSI = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DD_M2_COSI.jld", "grid")
grid_DDSTG_MM = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DDSTG_MM.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_DD_M2_COSI.y, grid_DD_M2_COSI.x, grid_DD_M2_COSI.value[:chisqr_cut] .- grid_DD_M2_COSI.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)

cs_DD_GAMMA_PBDOT = ax.contour(grid_DD_GAMMA_PBDOT.y, grid_DD_GAMMA_PBDOT.x, grid_DD_GAMMA_PBDOT.value[:chisqr_cut] .- grid_DD_GAMMA_PBDOT.params[:chisqr_min], levels=[lvl_68CL], colors="red")
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")

inds = filter(i -> (grid_DDSTG_MM.value[:chisqr_cut][i] .- grid_DDSTG_MM.params[:chisqr_min] <= lvl_68CL), 1:length(grid_DDSTG_MM.value[:chisqr]))
pclm2 = ax.plot(grid_DDSTG_MM.value[:M2][inds], sqrt.(1 .- grid_DDSTG_MM.value[:SINI][inds].^2), color="orange", ".")

ax.set_ylabel(get_label(tf.test.param1.name), size=16)
ax.set_xlabel(get_label(tf.test.param2.name), size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

#-------------------------------------------------------------------------------------\


grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_J1141-6545_MPA1_nice_90CL.jld", "grid")
grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_DDFWHE.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_DDSTG.y, grid_DDSTG.x, grid_DDSTG.value[:chisqr_cut] .- grid_DDSTG.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[0.0, 0.5, 1.0, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1], colors="black")
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[1.0, 1.5, 1.8], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

cs_DDSTG = ax.contour(grid_DDSTG.y, grid_DDSTG.x, grid_DDSTG.value[:chisqr] .- grid_DDSTG.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{DDSTG}", color=cm(3))

ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)

cs_DDFWHE = ax.contour(grid_DDFWHE.y, grid_DDFWHE.x, grid_DDFWHE.value[:chisqr] .- grid_DDFWHE.params[:chisqr_min], levels=levels=[lvl_90CL], colors=[cm(2)], linestyles=["-."])
plot([], [], label=L"\mathrm{PK\ method}", color=cm(2), ls="-.", zorder=4)

#title(L"\mathrm{PK\ method, 23-26}")
ax.set_ylabel(get_label(pkf.test.param1.name), size=16)
ax.set_xlabel(get_label(pkf.test.param2.name), size=16)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J1141-6545_DDSTG_chisqr_comp.pdf", dpi=400)
#-------------------------------------------------------------------------------------\

grid_DDSTG_J1141_6545 = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_J1141-6545_MPA1_nice_90CL.jld", "grid")
grid_PK_Double_Pulsar = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_Double_Pulsar.jld", "grid")
grid_PK_J2222_0137 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J2222−0137.jld", "grid")
grid_PK_J1738_0333 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J1738+0333.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()

cs_DDSTG_J1141_6545 = ax.contour(grid_DDSTG_J1141_6545.y, grid_DDSTG_J1141_6545.x, grid_DDSTG_J1141_6545.value[:chisqr_cut] .- grid_DDSTG_J1141_6545.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1952+2630", color=cm(3))


cs_Double_Pulsar = ax.contour(grid_PK_Double_Pulsar.y, grid_PK_Double_Pulsar.x, grid_PK_Double_Pulsar.value[:chisqr] .- grid_PK_Double_Pulsar.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(2)], linestyles=[(0, (3, 2, 1, 2, 1, 2))])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}0737-3039\mathrm{A}}", color=cm(2),ls=(0, (3, 2, 1, 2, 1, 2)))

cs_J2222_0137 = ax.contour(grid_PK_J2222_0137.y, grid_PK_J2222_0137.x, grid_PK_J2222_0137.value[:chisqr] .- grid_PK_J2222_0137.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(4)], linestyles=["-."])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}2222-0137", color=cm(4), ls="-.")

cs_J1738_0333 = ax.contour(grid_PK_J1738_0333.y, grid_PK_J1738_0333.x, grid_PK_J1738_0333.value[:chisqr] .- grid_PK_J1738_0333.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(0)], linestyles=[":"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1738-0333", color=cm(0), ls=":")


ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J1141-6545_pulsars_comp.pdf", dpi=400)