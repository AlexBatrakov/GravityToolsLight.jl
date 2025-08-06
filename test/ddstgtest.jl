using Revise
using GravityTools
using PyPlot
using JLD
pygui(true)

#-------------------------------------------------------------------------------------

sets = DDSTGTestSettings(
    par_file_init = "J2222-0137_T1_DDSTG_DMX.par",
    tim_file = "J2222-0137_T1.tim",
    add_flag = "-j -c",
    fit_XPBDOT = true)

ddstgt = DDSTGTest(
    psrname = "J2222-0137",
    eosname = "MPA1",
    N_refinement = 6,
    CL = 0.90,
    sets = sets)

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
ddstgt.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
calculate!(ddstgt)

#-------------------------------------------------------------------------------------

sets = DDSTGTestSettings(
    par_file_init = "J1141-6545_T1_DDSTG.par",
    tim_file = "J1141-6545_T1.tim",
    add_flag = "-c",
    fit_XPBDOT = true)

ddstgt = DDSTGTest(
    psrname = "J1141-6545",
    eosname = "MPA1",
    N_refinement = 6,
    CL = 0.90,
    sets = sets)

log10alpha0_grid = collect(LinRange(-4,-1,9))
beta0_grid = collect(LinRange(-6.0,+6.0,9))
grid_size_counter(grid_init=9,ref_level=6)
ddstgt.grid = SimpleGrid(Dict(), log10alpha0_grid, beta0_grid)
calculate!(ddstgt)

#-------------------------------------------------------------------------------------


cut_ddstg_grid!(ddstgt.grid)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(ddstgt.grid.y, ddstgt.grid.x, ddstgt.grid.value[:chisqr] .- ddstgt.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_5σ), rasterized=true)
cs = ax.contour(ddstgt.grid.y, ddstgt.grid.x, ddstgt.grid.value[:chisqr] .- ddstgt.grid.params[:chisqr_min], levels=[4.60517], colors="red")
#clabel(cs, cs.levels, fmt=Dict(lvl_2σ => "2σ"))
#cs.collections[1].set_label("$(ct.name)")
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}_{\mathrm{DDSTG}}$")
ax.axhline(y=log10(3.4e-3), color="gray", label="Cassini", ls="--", zorder=5)
plot([], [], label="90% CL", color="red")
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}|\alpha_0|", size=16)
legend(fontsize=12)
ylim(-4.0,-1.0)
tight_layout()
#savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_DDSTG_chisqr_map.pdf", dpi=600)


grid_ddstg = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves_old/grid_J2222-0137_contour_90CL.jld", "grid")
grid_ddstg = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")
grid_ddstg = load("/Users/abatrakov/Documents/PhD_work/projects/J1141-6545/saves/grid_J1141-6545_MPA1_nice_90CL.jld", "grid")
cut_ddstg_grid!(grid_ddstg)
rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
fig, ax = subplots()
pclm = ax.pcolormesh(grid_ddstg.y, grid_ddstg.x, grid_ddstg.value[:chisqr_cut] .- grid_ddstg.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_3σ), rasterized=true)
cs = ax.contour(grid_ddstg.y, grid_ddstg.x, grid_ddstg.value[:chisqr_cut] .- grid_ddstg.params[:chisqr_min], levels=[4.60517], colors="red")
#clabel(cs, cs.levels, fmt=Dict(lvl_2σ => "2σ"))
#cs.collections[1].set_label("$(ct.name)")
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}_{\mathrm{DDSTG}}$")
ax.axhline(y=log10(3.4e-3), color="gray", label="Cassini", ls="--", zorder=5)
plot([], [], label="DDSTG", color="red")
ax.set_xlabel(L"\beta_0", size=16)
ax.set_ylabel(L"\log_{10}|\alpha_0|", size=16)
legend(fontsize=12)
ylim(-4.0,-1.0)
tight_layout()
#savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/J2222-0137_DDSTG_chisqr_map.pdf", dpi=600)