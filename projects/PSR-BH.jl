using Revise
using GravityTools
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
pygui(true)

#-------------------------------------------------------------------------------------

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

toa_uncertainty = "1mus"

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
plot([], [], label=L"P_\mathrm{b} = 3 \mathrm{\ d}", color=cm(4), ls=":")

cs_1d = ax.contour(grid_1d.y, grid_1d.x, grid_1d.value[:chisqr_cut] .- grid_1d.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(6)], linestyles=[(0, (3, 2, 1, 2, 1, 2))])
#clabel(cs2, cs2.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"P_\mathrm{b} = 1 \mathrm{\ d}", color=cm(6), ls=(0, (3, 2, 1, 2, 1, 2)))

cs_8h = ax.contour(grid_8h.y, grid_8h.x, grid_8h.value[:chisqr_cut] .- grid_8h.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(0)])
#clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"P_\mathrm{b} = 8 \mathrm{\ h}", color=cm(0))



ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)

legend(fontsize=11)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$", size=14)
#title(L"0.1 \mu s")
#ax.invert_xaxis()
tight_layout()

savefig("/Users/abatrakov/Documents/PhD_work/projects/plots/PSRBH_$(toa_uncertainty)_chisqr_comp.pdf", dpi=400, bbox_inches="tight")
