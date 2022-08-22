using Revise
using GravityTools
using PyPlot
using JLD
using ColorSchemes
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
    alpha0 = -0.0,
    beta0 =  -0.0,
#    param1 = (name = "M2", min = 0.0, max = 2.0, N = 9),
    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
#    param2 = (name = "DTHETA", min = -500.0, max = 500.0, N = 9),
#    param1 = (name = "OMDOT", min = 5.3098, max = 5.3115, N = 9),
#    param1 = (name = "COSI", min = 0.0, max = 1.0, N = 9),
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
    N_refinement = 1,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = false
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
#    par_file_init = "J1952+2630_DDSTG_12jan21_new.par",
    par_file_init = "J1952+2630_DDSTG.par",
    tim_file = "TOAs_pdev_puppi_fast_T2_23-26_efac_gauss.tim",
#    tim_file = "TOAs_pdev_puppi_fast_T2_23-39_gauss.tim",
#    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
#    add_flag = "-c -j -ni npulses_23-26_efac.dat",
#    add_flag = "-c -j -ni npulses_23-39.dat",
#    add_flag = "-c -j",
    fit_XPBDOT = false,
    nits_first_step = 10,
    gain_fisrt_step = 1.0,
    params_first_step = [("m2", "PK", 0), ("mtot", "PK", 0), ("Pb", "PK", 1), ("e0", "PK", 1), ("x0", "PK", 1), ("omega0", "not changed", 0)],
    nits_second_step = 30,
    gain_second_step = 0.05,
    params_second_step = [("m2", "PK", 1), ("mtot", "PK", 1), ("Pb", "not changed", 1), ("e0", "not changed", 1), ("x0", "not changed", 1), ("omega0", "not changed", 1)]
    )

gsets = GridSetttings(
    N_refinement = 1,
    CL = [0.68],
    refinement_type = "contour",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = obs_params.masses_init.m1
pf.bnsys.comp.mass = obs_params.masses_init.m2
interpolate_bnsys!(pf)

calculate!(tf, pf, obs_params)

#-------------------------------------------------------------------------------------

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "J2222-0137_T1_DDSTG_DMX.par",
    tim_file = "J2222-0137_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = false,
    nits = 5,
    nits_conv = 0,
    gain_conv = 1.0,
    all_conv = false
    )

gsets = GridSetttings(
    N_refinement = 6,
    CL = 0.90,
    plot_type = "nice",
    chisqr_max = lvl_3σ,
    chisqr_delta = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)


#-------------------------------------------------------------------------------------

period = "1d"

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
    N_refinement = 6,
    CL = 0.90,
    plot_type = "nice",
    chisqr_max = lvl_3σ,
    chisqr_delta = 1.0,
    gr_in_chisqr = true
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
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
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

function get_label(name)
    if name == "PBDOT"
        return L"\dot{P}_\mathrm{b}\, (10^{-12})"
    elseif name == "M2"
        return L"m_{\mathrm{c}}\,(M_\odot)"
    elseif name == "MTOT"
        return L"m_{\mathrm{tot}}\,(M_\odot)"
    elseif name == "GAMMA"
        return L"\gamma"
    elseif name == "XDOT"
        return L"\dot{x}\,(10^{-12} s/s)"
    elseif name == "OMDOT"
        return L"\dot{\omega}"
    elseif name == "COSI"
        return L"\cos i"
    elseif name == "DTHETA"
        return L"\delta_\theta"
    elseif name == "log10alpha0"
        return L"\log_{10}|\alpha_0|"
    elseif name == "beta0"
        return L"\beta_0"
    elseif name == "H3"
        return L"h_3"
    elseif name == "VARSIGMA"
        return L"\varsigma"
    end
end

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
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=tf.gsets.delta_chisqr_max), rasterized=true)
cs = ax.contour(tf.grid.y, tf.grid.x, round.(tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], digits=1), levels=tf.gsets.contours, colors="red")
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")

cut_ddstg_grid!(tf.grid)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)
cs = ax.contour(tf.grid.y, tf.grid.x, round.(tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], digits=1), levels=[lvl_68CL], colors="red")
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")


#cs2 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:SINI], levels=[0.25, 0.50, 0.75,0.90, 0.95, 1.0], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=L"\sin i", color="black")
#cs2 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:XDOT], levels=[0.0, 0.15, 0.2, 0.25, 0.4, 0.5, 1.0], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=get_label("XDOT"), color="black")
#cs3 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:M1_GR], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="black")
#clabel(cs3, cs3.levels)
#plot([], [], label=L"m_{\mathrm{p,GR}}\,(M_\odot)", color="black")

cs4 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], levels=[1.0, lvl_68CL, lvl_95CL], colors="pink")
clabel(cs4, cs4.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"PK", color="pink")

cs_m1 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:MTOT] - tf.grid.value[:M2], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="black")
clabel(cs_m1, cs_m1.levels)
plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

cs_m2 = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:M2], levels=[0.0, 0.5, 1.0, 1.5, 2.0], colors="green")
clabel(cs_m2, cs_m2.levels)
plot([], [], label=L"m_{\mathrm{c}}\,(M_\odot)", color="green")

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
cs = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], levels=[lvl_68CL, lvl_90CL], colors="red")
clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")


#-------------------------------------------------------------------------------------

function get_smooth_contour(grid, lvl)
    cl = levels(contours(grid.x, grid.y, grid.value[:chisqr_cut] .- grid.params[:chisqr_min], [lvl]))[1]
    curves = []
    for line in lines(cl)
        Y, X = coordinates(line)
        piece = hcat(X,Y)
        piece_smooth = smooth_curve(piece, 50)
        curves = push!(curves, piece_smooth)
    end
    return curves
end

function plot_contour(ax, grid, lvl; color, label)
    smooth_contours = get_smooth_contour(grid, lvl)
    for (i, smooth_contour) in enumerate(smooth_contours)
        if i == 1
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color, label=label)
        else
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color)
        end
    end
end

plot_contour(ax, grid_obs,  lvl_90CL, color=cm(2), label="observations")