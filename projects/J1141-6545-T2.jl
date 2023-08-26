using Revise
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

using Distributed
addprocs(8)

@everywhere using GravityToolsLight

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# ln -s /Users/abatrakov/Documents/Work/PhD/computed_grids_fine data_ddstg
# cd /Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/
# cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/fake_with_gauss")
cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper")

test = GeneralTest(
    psrname = "J1141-6545",
    eosname = "BSk22",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 3),
    param2 = (name = "beta0", min = -6.0, max = 10.0, N = 3)
    )

# tsets = TempoSettings(
#     par_file_init = "J1141-6545_T1_DDSTG.par",
#     tim_file = "J1141-6545_T1.tim",
#     add_flag = "-nobs 22000 -newpar",
#     fit_XPBDOT = false,
#     nits_first_step = 3,
#     gain_fisrt_step = 1.0
#     )


# tempo2 -f J1141-6545_until_2018_DDSTG.par J1141-6545_until_2018.tim -nobs 22000 -newpar

tsets = TempoSettings(
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    add_flag = "-nobs 22000 -newpar",
    fit_XPBDOT = false,
    nits_first_step = 3,
    gain_fisrt_step = 1.0
    )

gsets = GridSetttings(
    N_refinement = 4,
    CL = [0.90],
    refinement_type = "contour",
    delta_chisqr_max = 10,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

# theory = DEF(0.0, 0.0)
# eosname = :MPA1
# bnsys = BinarySystem(:NS, :WD)
# sets = Settings("/Users/abatrakov/Documents/Science/Software/tempo2/T2runtime/ddstg_data")
# pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
# read_grid!(pf)

#obs_params = obs_params_dataset["J2222−0137_sim_DD"]
#pf.bnsys.K_params = obs_params.K
#pf.bnsys.psr.mass = obs_params.masses_init.m1
#pf.bnsys.comp.mass = obs_params.masses_init.m2
#interpolate_bnsys!(pf)

#calculate!(tf)

calculate_t2!(tf)
tf.grid.status

#calculate_t2!(tf, add_refinement=1)
#tf.grid.status



fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=tf.gsets.delta_chisqr_max), rasterized=true)
cs = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], levels=tf.gsets.contours, colors=["red"])

#-------------------------------------------------------------------------------------

grid_DDSTG_J1141_6545 = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/grid_DDSTG_J1141_6545_BSk22.jld", "grid")
grid_PK_J1141_6545 =    load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDFWHE_J1141-6545_BSk22.jld", "grid")
grid_PK_Triple_System = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Triple_System_BSk22.jld", "grid")
grid_PK_Double_Pulsar = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Double_Pulsar_BSk22.jld", "grid")
grid_PK_J1738_0333 =   load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_J1738+0333_Guo_BSk22.jld", "grid")

#grid_PK_J1141_6545_BHAT = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_J1141-6545_BHAT_68.jld", "grid")

# grid_PK_J2222_0137 = load("/Users/abatrakov/Documents/Work/PhD/projects/J1952+2630/saves/grid_PK_J2222−0137.jld", "grid")


rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()

pclm = ax.pcolormesh(grid_DDSTG_J1141_6545.y, grid_DDSTG_J1141_6545.x, grid_DDSTG_J1141_6545.value[:chisqr] .- grid_DDSTG_J1141_6545.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=6.0), rasterized=true)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")

# cs = ax.contour(grid_DDSTG_J1141_6545.y, grid_DDSTG_J1141_6545.x, grid_DDSTG_J1141_6545.value[:chisqr_cut] .- grid_DDSTG_J1141_6545.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(3)])
# clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
# plot([], [], label=L"\mathrm{J}1141-6545", color=cm(3))

plot_contour(ax, grid_DDSTG_J1141_6545,  lvl_90CL, color=cm(3), label=L"\mathrm{J}1141-6545\ \mathrm{DDSTG}", N_smooth=10)

cs_PK_J1141_6545 = ax.contour(grid_PK_J1141_6545.y, grid_PK_J1141_6545.x, grid_PK_J1141_6545.value[:chisqr] .- grid_PK_J1141_6545.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(0)], linestyles=[(0, (3, 2, 1, 2, 1, 2))])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1141-6545\ \mathrm{PK}", color=cm(0),ls=(0, (3, 2, 1, 2, 1, 2)))

# cs_PK_J1141_6545_BHAT = ax.contour(grid_PK_J1141_6545_BHAT.y, grid_PK_J1141_6545_BHAT.x, grid_PK_J1141_6545_BHAT.value[:chisqr] .- grid_PK_J1141_6545_BHAT.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(5)], linestyles=["-"])
# #clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
# plot([], [], label=L"\mathrm{J}1141-6545\ \mathrm{Bhat}", color=cm(5),ls="-")

cs_Double_Pulsar = ax.contour(grid_PK_Double_Pulsar.y, grid_PK_Double_Pulsar.x, grid_PK_Double_Pulsar.value[:chisqr] .- grid_PK_Double_Pulsar.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(4)], linestyles=["-."])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}0737-3039\mathrm{A}}", color=cm(4), ls="-.")

# cs_J2222_0137 = ax.contour(grid_PK_J2222_0137.y, grid_PK_J2222_0137.x, grid_PK_J2222_0137.value[:chisqr] .- grid_PK_J2222_0137.params[:chisqr_min], levels=[lvl_68CL], colors=[cm(4)], linestyles=["-."])
# #clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
# plot([], [], label=L"\mathrm{J}2222-0137", color=cm(4), ls="-.")


cs_J1738_0333 = ax.contour(grid_PK_J1738_0333.y, grid_PK_J1738_0333.x, grid_PK_J1738_0333.value[:chisqr] .- grid_PK_J1738_0333.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(2)], linestyles=[":"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1738+0333", color=cm(2), ls=":")

cs_J1738_0333 = ax.contour(grid_PK_Triple_System.y, grid_PK_Triple_System.x, grid_PK_Triple_System.value[:chisqr] .- grid_PK_Triple_System.params[:chisqr_min], levels=[lvl_90CL], colors=["red"], linestyles=["-"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{Triple System}", color="red", ls="-")


ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
#ax.invert_xaxis()
xlim(left=-6.0, right=+6.0)
ylim(top=-1.5, bottom=-4.0)
tight_layout()

savefig("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/plot_BSk22.pdf")