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

function plot_contour(ax, grid, lvl; color, label, N_smooth=20, ls="-", linewidth=1.5, zorder=nothing)
    smooth_contours = get_smooth_contour(grid, lvl, N_smooth)
    for (i, smooth_contour) in enumerate(smooth_contours)
        if i == 1
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color, label=label, ls=ls, linewidth=linewidth, zorder=zorder)
        else
            ax.plot(smooth_contour[:,1],smooth_contour[:,2], color=color, ls=ls, linewidth=linewidth, zorder=zorder)
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

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J2222−0137_Guo_DDK"]
obs_params = obs_params_dataset["J2222−0137_Guo_ELL1H+"]
obs_params = obs_params_dataset["J2222−0137_DDFWHE"]
obs_params = obs_params_dataset["J2222−0137_DDK"]
obs_params = obs_params_dataset["J2222−0137_DDK_XDOT"]
obs_params = obs_params_dataset["J2222−0137_DD"]
obs_params = obs_params_dataset["J2222−0137_sim_DD"]
obs_params = obs_params_dataset["J2222−0137_DDFWHE_bestXPBDOT"]

pf.bnsys.K_params = obs_params.K

pf.theory.alpha0 = -0.0
pf.theory.beta0  = 0.0
interpolate_mgrid!(pf)

pf.bnsys.psr.mass = 1.81
pf.bnsys.comp.mass = 1.31
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

optimize_PK_method(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J2222-0137",
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
    CL = [0.90],
    refinement_type = "nice",
    delta_chisqr_max = 25,
    delta_chisqr_diff = 0.2,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

#-------------------------------------------------------------------------------------

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings(ENV["TEMPO"] * "/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J2222−0137_DDFWHE"]
obs_params = obs_params_dataset["J2222−0137_DDK"]
obs_params = obs_params_dataset["J2222−0137_Guo_ELL1H+"]
obs_params = obs_params_dataset["J2222−0137_DD"]
#obs_params = obs_params_dataset["J2222−0137_DDK_XDOT"]
#obs_params = obs_params_dataset["J2222−0137_sim_DD"]

pf.bnsys.K_params = obs_params.K

pf.theory.alpha0 = -0.0
pf.theory.beta0  = 0.0
interpolate_mgrid!(pf)

pf.bnsys.psr.mass = 1.81
pf.bnsys.comp.mass = 1.31
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
#    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
#    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    alpha0 = -0.0001,
    beta0 =  -4.5,
#    param2 = (name = "PBDOT", min = -0.00, max = -0.03, N = 9),
#    param1 = (name = "OMDOT", min = 0.0963854 - 3*0.0004908, max = 0.0963854 +3*0.0004908, N = 9),
#    param2 = (name = "Pbdot", min = -0.32e-12, max = -0.45e-12, N = 9),
#    param1 = (name = "gamma", min = 0.0004, max = 0.0011, N = 9),
    param1 = (name = "m2", min = 0.0, max = 3.0, N = 9),
    param2 = (name = "m1", min = 0.0, max = 3.0, N = 9),
#    param2 = (name = "s", min = 0.1, max = 1.0, N = 9)
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

test = GeneralTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
#    par_file_init = "J2222-0137_T1_DDSTG_DMX.par",
#    tim_file = "J2222-0137_T1.tim",
#    add_flag = "-c -j -ni npulses.dat",
    par_file_init = "J2222-0137_DDSTG.par",
    tim_file = "epta_mk_fast.tim",
    add_flag = "-c -ni npulses.dat",
    fit_XPBDOT = false,
    nits_first_step = 5,
    gain_fisrt_step = 1.0
    )

gsets = GridSetttings(
    N_refinement = 5,
    CL = [0.90],
    refinement_type = "nice",
    delta_chisqr_max = 10,
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

obs_params = obs_params_dataset["J2222−0137_sim_DD"]
pf.bnsys.K_params = obs_params.K
pf.bnsys.psr.mass = obs_params.masses_init.m1
pf.bnsys.comp.mass = obs_params.masses_init.m2
interpolate_bnsys!(pf)

calculate!(tf, pf, obs_params)


#-------------------------------------------------------------------------------------

grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")

cut_ddstg_grid!(grid_DDSTG)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_DDSTG.y, grid_DDSTG.x, grid_DDSTG.value[:chisqr_cut] .- grid_DDSTG.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=25.0), rasterized=true)
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[0.0, 0.5, 1.0, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1], colors="black")
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[1.0, 1.5, 1.8], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

ax.axhline(y=log10(3.4e-3), color="grey", label=L"\mathrm{Cassini}", ls="--", zorder=5)

cs_DDSTG = ax.contour(grid_DDSTG.y, grid_DDSTG.x, grid_DDSTG.value[:chisqr] .- grid_DDSTG.params[:chisqr_min], levels=[lvl_90CL], colors=["red"])
#clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"90\mathrm{\%\ CL}", color="red")



#title(L"\mathrm{PK\ method, 23-26}")
ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
legend(fontsize=11)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_obs_chisqr_90CL.pdf", dpi=600, bbox_inches="tight")

#-------------------------------------------------------------------------------------

grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves_old/grid_J2222-0137_contour_90CL.jld", "grid")

cm = ColorMap(ColorSchemes.:OrRd_7.colors, 7)
fig2, ax2 = subplots()
pclm2 = ax2.pcolormesh(grid_DDSTG.y, grid_DDSTG.x, grid_DDSTG.ref_level, cmap=cm, rasterized=true)
cbar2 = colorbar(pclm2)
cbar2.set_label(L"\mathrm{level\ of\ refinement}", size=16)
#ax2.set_title("$(psr_name) TOAs; $EOS eos", size=16)
ax2.set_xlabel(L"\beta_0", size=16)
ax2.set_ylabel(L"\log_{10}|\alpha_0|", size=16)
tight_layout()
savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_ref_level.pdf", dpi=600, bbox_inches="tight")

#-------------------------------------------------------------------------------------

grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL_fixedXPBDOT.jld", "grid")
grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")
grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_DDFWHE.jld", "grid")
#grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_Guo_DDK.jld", "grid")
#grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_DDK_XDOT.jld", "grid")
grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_DD.jld", "grid")
grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_DD_bestXPBDOT.jld", "grid")
grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_DDFWHE_bestXPBDOT.jld", "grid")
#grid_DDSTG = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_sim.jld", "grid")
#grid_DDFWHE = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_sim_DD.jld", "grid")

inds = filter(i -> (grid_DDSTG.value[:ALPHA0][i] >= -1e0), 1:length(grid_DDSTG.value[:chisqr]))


fig, ax = subplots()
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
pclm = ax.scatter(grid_DDSTG.value[:chisqr][inds] .- grid_DDSTG.params[:chisqr_gr], grid_DDFWHE.value[:chisqr][inds] .- grid_DDFWHE.params[:chisqr_gr] .- (grid_DDSTG.value[:chisqr][inds] .- grid_DDSTG.params[:chisqr_gr]), alpha=0.5, linewidth=0.05, c=grid_DDSTG.value[:BETA0][inds], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"(\chi^2 - \chi^2_\mathrm{GR})_\mathrm{DDSTG}", size=16)
ax.set_ylabel(L"(\chi^2 - \chi^2_\mathrm{GR})_\mathrm{PK} - (\chi^2 - \chi^2_\mathrm{GR})_\mathrm{DDSTG}", size=16)
axvline(x=4.6, color="red", label=L"90\mathrm{\%\ CL}")
legend(fontsize=12)
xlim(-0.5,10.0)
ylim(-1.0,1.0)
tight_layout()
#savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_DDSTG_PK_comp.pdf", dpi=600)

fig, ax = subplots()
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
pclm = ax.scatter(grid_DDSTG.value[:chisqr][inds] .- grid_DDSTG.params[:chisqr_min], grid_DDFWHE.value[:chisqr][inds] .- grid_DDFWHE.params[:chisqr_min] .- (grid_DDSTG.value[:chisqr][inds] .- grid_DDSTG.params[:chisqr_min]), alpha=0.5, linewidth=0.05, c=grid_DDSTG.value[:BETA0][inds], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta\chi^2_\mathrm{DDSTG}", size=16)
ax.set_ylabel(L"\Delta\chi^2_\mathrm{PK} - \Delta\chi^2_\mathrm{DDSTG}", size=16)
axvline(x=4.6, color="red", label=L"90\mathrm{\%\ CL}")
legend(fontsize=12)
xlim(-0.5,10.0)
ylim(-3.0,0.5)
tight_layout()
savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_DDSTG_PK_DD_chisqr_comp.pdf", dpi=600, bbox_inches="tight")

fig, ax = subplots()
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
pclm = ax.scatter(grid_DDSTG.value[:chisqr][inds] .- grid_DDSTG.params[:chisqr_min], grid_DDFWHE.value[:chisqr][inds] .- grid_DDFWHE.params[:chisqr_min], alpha=0.5, linewidth=0.05, c=grid_DDSTG.value[:BETA0][inds], norm = matplotlib.colors.Normalize(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
cbar.set_label(L"$\beta_0$", size=16)
ax.set_xlabel(L"\Delta\chi^2_\mathrm{DDSTG}", size=16)
ax.set_ylabel(L"\Delta\chi^2_\mathrm{PK}", size=16)
axvline(x=4.6, color="red", label=L"90\mathrm{\%\ CL}")
axhline(y=4.6, color="red")
legend(fontsize=12)
xlim(-0.5,10.0)
ylim(-0.5,10.0)
tight_layout()


fig, ax = subplots()
pclm = ax.scatter(grid_ddstg.value[:chisqr][inds] .- grid_ddstg.params[:chisqr_min], ct.grid.value[:chisqr][inds] .- ct.grid.params[:chisqr_min] .- (grid_ddstg.value[:chisqr][inds] .- grid_ddstg.params[:chisqr_min]), alpha=0.3, linewidth=0.05, c=-grid_ddstg.value[:ALPHA0][inds], norm = matplotlib.colors.LogNorm(), rasterized=true)
cbar = colorbar(pclm)
cbar.set_alpha(1)
cbar.draw_all()
cbar.set_label(L"$|\alpha_0|$", size=16)
ax.set_xlabel(L"\Delta \chi^2_{\mathrm{DDSTG}}", size=16)
ax.set_ylabel(L"\Delta \chi^2_{\mathrm{PK}} - \Delta \chi^2_{\mathrm{DDSTG}}", size=16)
axvline(x=4.6, color="red", label="90% CL")
legend(fontsize=12)
xlim(0.0,10.0)
ylim(-2.0,0.25)
tight_layout()
#savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_PK_DDSTG_comp_chisqr_alpha0_cut.pdf", dpi=600)
#-------------------------------------------------------------------------------------

grid_obs = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_contour_90CL.jld", "grid")
grid_all = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_all_nice_90CL.jld", "grid")
grid_fast = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_fast_contour_90CL.jld", "grid")
grid_mk = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_mk_contour_90CL.jld", "grid")
grid_epta = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_epta_contour_90CL.jld", "grid")

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
axhline(y=log10(3.4e-3), color="gray", label=L"\mathrm{Cassini}", ls="--")
plot_contour(ax, grid_obs,  lvl_90CL, color=cm(2), label=L"\mathrm{Real\ data}", N_smooth=50, ls="-.")
plot_contour(ax, grid_epta, lvl_90CL, color=cm(6), label=L"\mathrm{3ERT}", N_smooth=50, ls=(0, (3, 1, 1, 1)))
plot_contour(ax, grid_mk,   lvl_90CL, color=cm(0), label=L"\mathrm{MK}", N_smooth=50, ls=(0, (3, 2, 1, 2, 1, 2)))
plot_contour(ax, grid_fast, lvl_90CL, color=cm(4), label=L"\mathrm{FAST}", N_smooth=80, ls=":", linewidth=2.0, zorder=5)
plot_contour(ax, grid_all,  lvl_90CL, color=cm(1), label=L"\mathrm{FAST+MK+3ERT}", N_smooth=80)
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}(|\alpha_0|)", size=16)
ylim(-4.0,-1.0)
xlim(-6.0,6.0)
legend(fontsize=11, loc=1)
tight_layout()
savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_sim_chisqr_comp_90CL.pdf", dpi=600)

#-------------------------------------------------------------------------------------

grid_obs = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_contour_90CL.jld", "grid")
grid_sim = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_all_nice_90CL.jld", "grid")

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_sim.y, grid_sim.x, grid_sim.value[:chisqr_cut] .- grid_sim.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=25), rasterized=true)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
ax.axhline(y=log10(3.4e-3), color="gray", label=L"\mathrm{Cassini}", zorder=5, ls="--")
plot_contour(ax, grid_obs,  lvl_90CL, color=cm(2), label=L"\mathrm{Real\ data}", N_smooth=50, ls="-.")
plot_contour(ax, grid_sim, lvl_90CL, color=cm(0), label=L"\mathrm{FAST+MK+3ERT}", N_smooth=80)
#cs1 = ax.contour(grid_sim.y, grid_sim.x, grid_sim.value[:chisqr_cut] .- grid_sim.params[:chisqr_min], levels=[lvl_90CL], colors="black")
#cs1.collections[1].set_label("simulation")
#ax.set_title("$(psr_name) TOAs; $EOS eos")
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}(|\alpha_0|)", size=16)
legend(fontsize=11)
ylim(-4.0,-1.0)
tight_layout()
savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_sim_chisqr_90CL.pdf", dpi=600)

#-------------------------------------------------------------------------------------

grid_obs = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_contour_90CL.jld", "grid")
grid_all = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_all_nice_90CL.jld", "grid")
grid_fast = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_fast_contour_90CL.jld", "grid")
grid_mk = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_mk_contour_90CL.jld", "grid")
grid_epta = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/grid_J2222-0137_sim_sep_epta_contour_90CL.jld", "grid")

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_all.y, grid_all.x, grid_all.value[:chisqr_cut] .- grid_all.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=25), rasterized=true)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
ax.axhline(y=log10(3.4e-3), color="gray", label=L"\mathrm{Cassini}", zorder=5, ls="--")
plot_contour(ax, grid_obs,  lvl_90CL, color=cm(2), label=L"\mathrm{Real\ data}", N_smooth=50, ls="-.")
plot_contour(ax, grid_epta, lvl_90CL, color=cm(1), label=L"\mathrm{3ERT}", N_smooth=50, ls=(0, (3, 1, 1, 1)))
plot_contour(ax, grid_mk,   lvl_90CL, color=cm(4), label=L"\mathrm{MK}", N_smooth=50, ls=(0, (3, 2, 1, 2, 1, 2)))
plot_contour(ax, grid_fast, lvl_90CL, color=cm(0), label=L"\mathrm{FAST}", N_smooth=80, ls=":", linewidth=2.0, zorder=5)
plot_contour(ax, grid_all,  lvl_90CL, color=cm(6), label=L"\mathrm{FAST+MK+3ERT}", N_smooth=80)
#cs1 = ax.contour(grid_sim.y, grid_sim.x, grid_sim.value[:chisqr_cut] .- grid_sim.params[:chisqr_min], levels=[lvl_90CL], colors="black")
#cs1.collections[1].set_label("simulation")
#ax.set_title("$(psr_name) TOAs; $EOS eos")
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}|\alpha_0|", size=16)
legend(fontsize=10, loc=1)
ylim(-4.0,-1.0)
tight_layout()
savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_sim_chisqr_90CL_all.pdf", dpi=600,bbox_inches="tight")