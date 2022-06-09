using Revise
using GravityTools
using PyPlot
using JLD
pygui(true)


#test = GeneralTest(
#    psrname = "J1141-6545",
#    eosname = "MPA1",
#    alpha0  = -1e-4,
#    param1 = (name = "beta0", min = -6.0, max = 6.0, N = 9),
#    param2 = (name = "XPBDOT", min = -1e-13, max = 1e-13, N = 9)
#    )

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "J1141-6545_T1_DDSTG.par",
    tim_file = "J1141-6545_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = true
    )

gsets = GridSetttings(
    N_refinement = 1,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J1141-6545",
#    param1 = (name = "M2", min = 0.0, max = 2.0, N = 9),
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
#    param2 = (name = "DTHETA", min = -500.0, max = 500.0, N = 9),
#    param1 = (name = "OMDOT", min = 5.3098, max = 5.3115, N = 9),
#    param1 = (name = "COSI", min = 0.0, max = 1.0, N = 9),
    param1 = (name = "H3", min = 0.0, max = 0.000006, N = 9),
#    param1 = (name = "VARSIGMA", min = 0.0, max = 1.5, N = 9),
    param2 = (name = "XDOT", min = -1.0, max = 1.5, N = 9)
    )

tsets = TempoSettings(
#    par_file_init = "J1141-6545_T1_DD.par",
    par_file_init = "J1141-6545_T1_DDFWHE.par",
    tim_file = "J1141-6545_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = false,
    nits = 5,
    nits_conv = 10,
    gain_conv = 0.1,
    all_conv = true
    )

gsets = GridSetttings(
    N_refinement = 3,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J1952+2630",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "J1952+2630_DDSTG_12jan21_new.par",
    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false,
    nits = 0,
    nits_conv = 400,
    gain_conv = 0.03,
    all_conv = true
    )

gsets = GridSetttings(
    N_refinement = 0,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

#-------------------------------------------------------------------------------------

period = "3d"

test = GeneralTest(
    psrname = "B2127+11C$period",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "DDSTG_PSRBH$period.par",
    tim_file = "DDSTG_PSRBH$period.tim",
    add_flag = "-c",
    fit_XPBDOT = true,
    nits = 10,
    nits_conv = 0,
    gain_conv = 1.0,
    all_conv = false
    )

gsets = GridSetttings(
    N_refinement = 2,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J1952+2630",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "J1952+2630_DDSTG_12jan21_new.par",
    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false,
    nits = 0,
    nits_conv = 400,
    gain_conv = 0.03,
    all_conv = true
    )

gsets = GridSetttings(
    N_refinement = 0,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo-13.103_ddstg_arm_clock/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = 1.241308
pf.bnsys.comp.mass = 0.950913
interpolate_bnsys!(pf)

calculate!(tf, pf, obs_params)

#-------------------------------------------------------------------------------------
test = GeneralTest(
    psrname = "J1952+2630",
    eosname = "MPA1",
    alpha0 = -0.001,
    beta0 =  -6.0,
    param1 = (name = "MTOT", min = 0.5, max = 3.0, N = 9),
    param2 = (name = "M2", min = 0.5, max = 2.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "J1952+2630_DDSTG_12jan21_new.par",
    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
    fit_XPBDOT = false
    )

gsets = GridSetttings(
    N_refinement = 1,
    CL = 0.90,
    plot_type = "nice"
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

#-------------------------------------------------------------------------------------


function calculate_m1_gr!(tf)
    tf.grid.value[:SINI] = 2 .* tf.grid.value[:VARSIGMA] ./ (1 .+ tf.grid.value[:VARSIGMA] .^2)
    tf.grid.value[:M2] = tf.grid.value[:H3] ./ tf.grid.value[:VARSIGMA] .^3 * c^3 / G_CAV ./ M_sun
    tf.grid.value[:M1_GR] = @. ( tf.grid.value[:SINI] * tf.grid.value[:M2]*M_sun * (G_CAV*2*pi/tf.grid.value[:PB]/24/3600)^(1/3) / (c*2*pi/tf.grid.value[:PB]/24/3600*tf.grid.value[:A1]))^1.5 ./ M_sun - tf.grid.value[:M2]
end

calculate_m1_gr!(tf)
cut_ddstg_grid!(tf.grid)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_5σ), rasterized=true)
cs = ax.contour(tf.grid.y, tf.grid.x, round.(tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], digits=1), levels=[lvl_68CL, lvl_95CL], colors="red")
clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")
#cs2 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:SINI], levels=[0.25, 0.50, 0.75,0.90, 0.95, 1.0], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=L"\sin i", color="black")
#cs2 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:XDOT], levels=[0.0, 0.15, 0.2, 0.25, 0.4, 0.5, 1.0], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=get_label("XDOT"), color="black")
cs3 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:M1_GR], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="black")
clabel(cs3, cs3.levels)
plot([], [], label=L"m_{\mathrm{p,GR}}\,(M_\odot)", color="black")

cs4 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="green")
clabel(cs4, cs4.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"PK", color="green")

cs5 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:MTOT] - tf.grid.value[:M2], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="black")
clabel(cs5, cs5.levels)
plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

#clabel(cs, cs.levels, fmt=Dict(lvl_2σ => "2σ"))
#cs.collections[1].set_label("$(ct.name)")
ax.set_ylabel(get_label(tf.test.param1.name), size=16)
ax.set_xlabel(get_label(tf.test.param2.name), size=16)
ax.axhline(y=log10(3.4e-3), color="gray", label="Cassini", ls="--", zorder=5)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()


fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:MTOT] .- tf.grid.value[:M2], cmap="Blues_r", norm = matplotlib.colors.Normalize(), rasterized=true)
cs = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="red")
clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")

#-------------------------------------------------------------------------------------

tf_8h
tf_1d
tf_3d

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf_8h.grid.y, tf_8h.grid.x, tf_8h.grid.value[:chisqr] .- tf_8h.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_5σ), rasterized=true)

ax.axhline(y=log10(3.4e-3), color="gray", label=L"\mathrm{Cassini}", ls="--", zorder=5)

cs3 = ax.contour(tf_3d.grid.y, tf_3d.grid.x, tf_3d.grid.value[:chisqr] .- tf_3d.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="green")
clabel(cs3, cs3.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"3 \mathrm{\ days}", color="green")

cs2 = ax.contour(tf_1d.grid.y, tf_1d.grid.x, tf_1d.grid.value[:chisqr] .- tf_1d.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="blue")
clabel(cs2, cs2.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"1 \mathrm{\ day}", color="blue")

cs = ax.contour(tf_8h.grid.y, tf_8h.grid.x, tf_8h.grid.value[:chisqr] .- tf_8h.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_95CL], colors="red")
clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"8 \mathrm{\ hours}", color="red")


ax.set_ylabel(get_label(tf_8h.test.param1.name), size=16)
ax.set_xlabel(get_label(tf_8h.test.param2.name), size=16)

legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
title(L"0.1 \mu s")
#ax.invert_xaxis()
tight_layout()