using Revise
using GravityTools
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

function smooth_curve(curve, n)
    curve_smooth = copy(curve)
    n_half = div(n,2)
    N = length(curve[:,1])
    ext_curve_1 = vcat(2*curve[1,1] .- curve[n_half+1:-1:2,1], curve[:,1], 2*curve[N,1] .- curve[N-1:-1:N-1-n_half,1])
    ext_curve_2 = vcat(2*curve[1,2] .- curve[n_half+1:-1:2,2], curve[:,2], 2*curve[N,2] .- curve[N-1:-1:N-1-n_half,2])
    for i in 1:N
        inds = (i-n_half) + n_half : (i+n_half) + n_half
        curve_smooth[i,1] = mean(ext_curve_1[inds])
        curve_smooth[i,2] = mean(ext_curve_2[inds]) 
    end
    return curve_smooth
end

function get_smooth_contour(grid, lvl, N_smooth=20)
    cl = levels(contours(grid.x, grid.y, grid.value[:chisqr_cut] .- grid.params[:chisqr_min], [lvl]))[1]
    curves = []
    for line in lines(cl)
        Y, X = coordinates(line)
        if length(Y) - 1 - div(N_smooth,2) > 1
            piece = hcat(X,Y)
            piece_smooth = smooth_curve(piece, N_smooth)
            curves = push!(curves, piece_smooth)
        end
    end
    return curves
end

function plot_contour(ax, grid, lvl; color, label, N_smooth=20)
    smooth_contours = get_smooth_contour(grid, lvl, N_smooth)
    for (i, smooth_contour) in enumerate(smooth_contours)
        if i == 1
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color, label=label)
        else
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color)
        end
    end
end

function compare_contours(grid, lvl)
    rc("mathtext",fontset="cm")
    rc("font", family="serif", size=12)
    fig, ax = subplots()
    pclm = ax.pcolormesh(grid.y, grid.x, grid.value[:chisqr_cut] .- grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)
    cbar = colorbar(pclm)
    cbar.set_label(L"$\Delta\chi^{2}$")
    cs = ax.contour(grid.y, grid.x, grid.value[:chisqr_cut] .- grid.params[:chisqr_min], levels=[lvl], colors="red", label="orig")
    plot_contour(ax, grid,  lvl, color="black", label="smooth")
    axhline(y=log10(3.4e-3), color="gray", label="Cassini")
    ax.set_xlabel(L"\beta_0", size=16)
    ax.set_ylabel(L"\log_{10}(|\alpha_0|)", size=16)
    legend(fontsize=12)
    ylim(-4.0,-1.0)
end

function get_label(name)
    if name == "PBDOT"
        return L"\dot{P}_\mathrm{b}\, (10^{-12} s/s)"
    elseif name == "Pbdot"
        return L"\dot{P}_\mathrm{b}\, (s/s)"
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
#    tim_file = "TOAs_pdev_puppi_fast_T2_23-39_efac_gauss.tim",
#    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
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
    N_refinement = 5,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = obs_params.masses_init.m1
pf.bnsys.comp.mass = obs_params.masses_init.m2
interpolate_bnsys!(pf)

calculate!(tf, pf, obs_params)
save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_23_26_efac.jld", "grid", tf.grid)

#-------------------------------------------------------------------------------------

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
#    tim_file = "TOAs_pdev_puppi_fast_T2_23-26_efac_gauss.tim",
#    tim_file = "TOAs_pdev_puppi_fast_T2_23-39_efac_gauss.tim",
    tim_file = "TOAs_4o_2032_gauss.tim",
#    tim_file = "TOA_puppi_fast.tim",
    add_flag = "-c -j -ni npulses.dat",
#    add_flag = "-c -j",
    fit_XPBDOT = false,
    nits_first_step = 10,
    gain_fisrt_step = 1.0,
    params_first_step = [("m2", "PK", 0), ("mtot", "PK", 0), ("Pb", "PK", 1), ("e0", "PK", 1), ("x0", "PK", 1), ("omega0", "not changed", 0), ("XPBDOT", "0.0069", 0)],
    nits_second_step = 30,
    gain_second_step = 0.05,
    params_second_step = [("m2", "PK", 1), ("mtot", "PK", 1), ("Pb", "not changed", 1), ("e0", "not changed", 1), ("x0", "not changed", 1), ("omega0", "not changed", 1), ("XPBDOT", "0.0069", 0)]
    )

gsets = GridSetttings(
    N_refinement = 5,
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
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_1SXPBDOT"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = obs_params.masses_init.m1
pf.bnsys.comp.mass = obs_params.masses_init.m2
interpolate_bnsys!(pf)

calculate!(tf, pf, obs_params)

#save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_4o_32.jld", "grid", tf.grid)

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-39_efac"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_XPBDOT"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_1SXPBDOT"] 
obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_1.5M_1SXPBDOT"]
obs_params = obs_params_dataset["J1738+0333"]
#obs_params = obs_params_dataset["J2222−0137_Guo_DDK"]
#obs_params = obs_params_dataset["Triple System"]
#obs_params = obs_params_dataset["Double Pulsar"]
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
    eosname = "MPA1",
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
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-39_efac"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_XPBDOT"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_1SXPBDOT"] 
#obs_params = obs_params_dataset["J1738+0333"]
#obs_params = obs_params_dataset["J2222−0137_Guo_DDK"]
#obs_params = obs_params_dataset["Triple System"]
#obs_params = obs_params_dataset["Double Pulsar"]
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
    eosname = "MPA1",
    alpha0 = -0.00,
    beta0 =  -0.00,
#    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
#    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
    param1 = (name = "m2", min = 0.5, max = 1.5, N = 9),
    param2 = (name = "m1", min = 0.5, max = 2.0, N = 9),
    )

gsets = GridSetttings(
    N_refinement = 7,
    CL = [0.68],
    refinement_type = "massmass",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)
#-------------------------------------------------------------------------------------


save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_real_contour.jld", "grid", tf.grid)
save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_23_26.jld", "grid", tf.grid)
save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_23_39_efac.jld", "grid", tf.grid)

save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_Double_Pulsar.jld", "grid", pkf.grid)
save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J2222−0137.jld", "grid", pkf.grid)
save("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J1738+0333.jld", "grid", pkf.grid)

#-------------------------------------------------------------------------------------

grid_real = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_real_contour.jld", "grid")
grid_23_26 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_23_26_efac.jld", "grid")
grid_23_39 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_23_39_efac.jld", "grid")
grid_4o_32 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_4o_32_new.jld", "grid")
grid_4o_32_1SXPBDOT = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_4o_32_1SXPBDOT_new.jld", "grid")

cut_ddstg_grid!(grid_23_26)
cut_ddstg_grid!(grid_23_39)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_4o_32.y, grid_4o_32.x, grid_4o_32.value[:chisqr_cut] .- grid_4o_32.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)

#cs_real = ax.contour(grid_real.y, grid_real.x, grid_real.value[:chisqr_cut] .- grid_real.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(4)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"real", color=cm(4))

plot_contour(ax, grid_real,  lvl_68CL, color=cm(0), label=L"\mathrm{Real\ data}", N_smooth=30)

#cs_23_26 = ax.contour(grid_23_26.y, grid_23_26.x, grid_23_26.value[:chisqr_cut] .- grid_23_26.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(4)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"2023 - 2026", color=cm(4))

#cs_23_39 = ax.contour(grid_23_39.y, grid_23_39.x, grid_23_39.value[:chisqr_cut] .- grid_23_39.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"2023 - 2039", color=cm(3))

#cs_grid_4o_32_1SXPBDOT = ax.contour(grid_4o_32_1SXPBDOT.y, grid_4o_32_1SXPBDOT.x, grid_4o_32_1SXPBDOT.value[:chisqr_cut] .- grid_4o_32.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(2)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"\mathrm{Gal\ corr}", color=cm(2))

plot_contour(ax, grid_4o_32_1SXPBDOT,  lvl_68CL, color=cm(2), label=L"\mathrm{Gal\ corr.}", N_smooth=30)


#cs_grid_4o_32 = ax.contour(grid_4o_32.y, grid_4o_32.x, grid_4o_32.value[:chisqr_cut] .- grid_4o_32.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(1)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"\mathrm{No\ Gal\ corr}", color=cm(1))

plot_contour(ax, grid_4o_32,  lvl_68CL, color=cm(1), label=L"\mathrm{No\ Gal\ corr.}", N_smooth=30)

ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12, loc=1)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}_\mathrm{DDSTG}$", size=14)
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J1952+2630_DDSTG_chisqr_comp.pdf", dpi=400)

#-------------------------------------------------------------------------------------

calculate_OMDOT_GR(grid_23_39)

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_23_39.y, grid_23_39.x, abs.((grid_23_39.value[:OMDOT] .- grid_23_39.value[:OMDOT_GR]) ./ grid_23_39.value[:OMDOT_GR]), cmap="Blues_r", rasterized=true, norm=matplotlib.colors.LogNorm(vmin=1e-7))

cs_23_39 = ax.contour(grid_23_39.y, grid_23_39.x, grid_23_39.value[:chisqr_cut] .- grid_23_39.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"2023 - 2039", color=cm(3))

ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$|\dot{\omega} - \dot{\omega}_{GR}| / \dot{\omega}_{GR}$")
#ax.invert_xaxis()
tight_layout()


function calculate_OMDOT_GR(grid)
    Pb = grid.value[:PB] .* d
    n = 2*pi ./ Pb
	e = grid.value[:E]
	m = grid.value[:MTOT] .* M_sun
    grid.value[:OMDOT_GR] = @. 3*((G_CAV*m*n)/c^3)^(2/3)/(1 - e^2) * 360 * 365.25 * d / Pb
end

#-------------------------------------------------------------------------------------

fig, ax = subplots()
pclm = ax.scatter(grid_4o_32_1SXPBDOT.value[:chisqr] .- grid_4o_32_1SXPBDOT.params[:chisqr_min], abs.(grid_4o_32_1SXPBDOT.value[:ALPHAA] .- grid_4o_32_1SXPBDOT.value[:ALPHA0]), alpha=0.2, linewidth=0.05, c=grid_4o_32_1SXPBDOT.value[:BETA0], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
#ax.set_yscale("symlog", linthresh=1e-3)
#ax.set_yscale("log")
ylim(0.0, 0.002)
xlim(-0.05,5.0)
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta \chi^2_{\mathrm{DDSTG}}", size=16)
ax.set_ylabel(L"|\alpha_\mathrm{p} - \alpha_\mathrm{c}|", size=16)
axvline(x=2.28, color="red", label=L"68\mathrm{\%\ CL}")
legend(fontsize=12)
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J1952+2630_alphap-alphac.pdf", dpi=600)

#-------------------------------------------------------------------------------------

fig, ax = subplots()
pclm = ax.scatter(grid_23_39.value[:chisqr] .- grid_23_39.params[:chisqr_min], grid_23_39.value[:MTOT] - grid_23_39.value[:M2], alpha=0.1, linewidth=0.05, c=grid_23_39.value[:BETA0], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
ylim(1.0, 1.3)
xlim(-0.2,10.0)
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta \chi^2_{\mathrm{DDSTG}}", size=16)
ax.set_ylabel(L"m_{p} (M_\odot)", size=16)
legend(fontsize=12)
tight_layout()

#-------------------------------------------------------------------------------------

fig, ax = subplots()
pclm = ax.scatter(grid_23_39.value[:chisqr] .- grid_23_39.params[:chisqr_min], grid_23_39.value[:M2], alpha=0.1, linewidth=0.05, c=grid_23_39.value[:BETA0], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
ylim(0.8, 1.0)
xlim(-0.2,10.0)
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta \chi^2_{\mathrm{DDSTG}}", size=16)
ax.set_ylabel(L"m_{c} (M_\odot)", size=16)
legend(fontsize=12)
tight_layout()

#-------------------------------------------------------------------------------------

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_23_39.y, grid_23_39.x, grid_23_39.value[:PBDOT], cmap="Blues_r", rasterized=true)

cs_23_39 = ax.contour(grid_23_39.y, grid_23_39.x, grid_23_39.value[:chisqr_cut] .- grid_23_39.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"2023 - 2039", color=cm(3))

ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\dot{P}_b$")
#ax.invert_xaxis()
tight_layout()


#-------------------------------------------------------------------------------------

grid_4o_32_1SXPBDOT = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_4o_32_1SXPBDOT_new.jld", "grid")
grid_PK_Double_Pulsar = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_Double_Pulsar.jld", "grid")
grid_PK_J2222_0137 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J2222−0137.jld", "grid")
grid_PK_J1738_0333 = load("/Users/abatrakov/Documents/PhD_work/projects/J1952+2630/saves/grid_PK_J1738+0333.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()

#cs_4o_32_1SXPBDOT = ax.contour(grid_4o_32_1SXPBDOT.y, grid_4o_32_1SXPBDOT.x, grid_4o_32_1SXPBDOT.value[:chisqr_cut] .- grid_4o_32_1SXPBDOT.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(3)])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
#plot([], [], label=L"\mathrm{J}1952+2630", color=cm(3))

plot_contour(ax, grid_4o_32_1SXPBDOT,  lvl_68CL, color=cm(3), label=L"\mathrm{J}1952+2630", N_smooth=30)


cs_Double_Pulsar = ax.contour(grid_PK_Double_Pulsar.y, grid_PK_Double_Pulsar.x, grid_PK_Double_Pulsar.value[:chisqr] .- grid_PK_Double_Pulsar.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(2)], linestyles=[(0, (3, 2, 1, 2, 1, 2))])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}0737-3039\mathrm{A}}", color=cm(2),ls=(0, (3, 2, 1, 2, 1, 2)))

cs_J2222_0137 = ax.contour(grid_PK_J2222_0137.y, grid_PK_J2222_0137.x, grid_PK_J2222_0137.value[:chisqr] .- grid_PK_J2222_0137.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(4)], linestyles=["-."])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}2222-0137", color=cm(4), ls="-.")

cs_J1738_0333 = ax.contour(grid_PK_J1738_0333.y, grid_PK_J1738_0333.x, grid_PK_J1738_0333.value[:chisqr] .- grid_PK_J1738_0333.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(0)], linestyles=[":"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1738+0333", color=cm(0), ls=":")


ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J1952+2630_pulsars_comp.pdf", dpi=400)

