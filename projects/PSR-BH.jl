using Revise
using GravityTools
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

#-------------------------------------------------------------------------------------


toa_uncertainty = "1mus"
period = "1d"

for toa_uncertainty in ["1mus", "10mus"], period in ["8h", "1d", "3d"]

cd("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/PSRBH_simulation/work_data_ram/$(toa_uncertainty)/")

test = GeneralTest(
    psrname = "B2127+11C$period",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
    )

tsets = TempoSettings(
    par_file_init = "DDSTG_PSRBH$period.par",
    tim_file = "DDSTG_PSRBH$period.tim",
    add_flag = "-c -ni npulses.dat",
    fit_XPBDOT = true,
    nits_first_step = 10,
    gain_fisrt_step = 1.0
    )

gsets = GridSetttings(
    N_refinement = 5,
    CL = [0.90],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

tf = TempoFramework(test, tsets, gsets)

calculate!(tf)

save("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/saves/grid_PSRBH_$(toa_uncertainty)_$period.jld", "grid", tf.grid)

end


#-------------------------------------------------------------------------------------

save("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/saves/grid_PSRBH_$(toa_uncertainty)_$period.jld", "grid", tf.grid)

#-------------------------------------------------------------------------------------

toa_uncertainty = "10mus"

grid_8h = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/saves/grid_PSRBH_$(toa_uncertainty)_8h.jld", "grid")
grid_1d = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/saves/grid_PSRBH_$(toa_uncertainty)_1d.jld", "grid")
grid_3d = load("/Users/abatrakov/Documents/PhD_work/projects/Huanchen/saves/grid_PSRBH_$(toa_uncertainty)_3d.jld", "grid")
grid_J2222_0137 = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
#cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_8h.y, grid_8h.x, grid_8h.value[:chisqr_cut] .- grid_8h.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0), rasterized=true)

ax.axhline(y=log10(3.4e-3), color="gray", label=L"\mathrm{Cassini}", ls="--", zorder=5)

cs_J2222_0137 = ax.contour(grid_J2222_0137.y, grid_J2222_0137.x, grid_J2222_0137.value[:chisqr] .- grid_J2222_0137.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(2)], linestyles=["-."])
plot([], [], label=L"\mathrm{J2222-0137}", color=cm(2), ls="-.", zorder=4)

cs_3d = ax.contour(grid_3d.y, grid_3d.x, grid_3d.value[:chisqr_cut] .- grid_3d.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(4)], linestyles=[":"])
#clabel(cs3, cs3.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\dot{P}_\mathrm{b} = 3 \mathrm{\ days}", color=cm(4), ls=":")

cs_1d = ax.contour(grid_1d.y, grid_1d.x, grid_1d.value[:chisqr_cut] .- grid_1d.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(6)], linestyles=[(0, (3, 2, 1, 2, 1, 2))])
#clabel(cs2, cs2.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\dot{P}_\mathrm{b} = 1 \mathrm{\ day}", color=cm(6), ls=(0, (3, 2, 1, 2, 1, 2)))

cs_8h = ax.contour(grid_8h.y, grid_8h.x, grid_8h.value[:chisqr_cut] .- grid_8h.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(0)])
#clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\dot{P}_\mathrm{b} = 8 \mathrm{\ hours}", color=cm(0))



ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)

legend(fontsize=11)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$", size=14)
#title(L"0.1 \mu s")
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/PSRBH_$(toa_uncertainty)_chisqr_comp.pdf", dpi=400)
