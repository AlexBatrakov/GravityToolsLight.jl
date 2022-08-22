using Revise
using GravityTools
using PyPlot
using JLD
using ColorSchemes
pygui(true)

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
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
    N_refinement = 4,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J2222−0137_Guo_DDK"]
pf.bnsys.K_params = obs_params.K

pf.bnsys.psr.mass = 1.81
pf.bnsys.comp.mass = 1.31
interpolate_bnsys!(pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

find_best_masses(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

pf.theory.alpha0 = -0.000100000
pf.theory.beta0  = -4.4
interpolate_mgrid!(pf)

optimize_PK_method(obs_params, pf)
GravityTools.get_chisqr(obs_params, pf)
check_terms_in_chisqr(obs_params, pf)

test = GeneralTest(
    psrname = "J2222−0137",
    eosname = "MPA1",
    param1 = (name = "log10alpha0", min = -4.0, max = -1.0, N = 9),
    param2 = (name = "beta0", min = -6.0, max = 6.0, N = 9)
#    param2 = (name = "PBDOT", min = -0.32, max = -0.45, N = 9),
#    param1 = (name = "GAMMA", min = 0.0004, max = 0.0011, N = 9),
    )

gsets = GridSetttings(
    N_refinement = 4,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)

#-------------------------------------------------------------------------------------\
theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1952+2630_DDFWHE"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-26_efac"]
#obs_params = obs_params_dataset["J1952+2630_DDFWHE_23-39"]
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
    N_refinement = 4,
    CL = [0.68],
    refinement_type = "nice",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = true
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)
#-------------------------------------------------------------------------------------

theory = DEF(0.0, 0.0)
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/Software/tempo_ram/data_ddstg")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)
read_grid!(pf)

obs_params = obs_params_dataset["J1141-6545_DD"]
obs_params = obs_params_dataset["J1952+2630_DDFWHE_4o_32_XPBDOT"]

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
    alpha0 = -0.0,
    beta0  = -0.0,
    param1 = (name = "Pbdot", min = -0.3863725e-12 * 0.90, max = -0.3863725e-12 * 1.1, N = 9),
#    param1 = (name = "m2", min = 1.016393 * 0.9, max = 1.016393 * 1.1, N = 9),
#    param1 = (name = "m2", min = 0.1, max = 2.0, N = 9),
#    param2 = (name = "m1", min = 1.273322 * 0.9, max = 1.273322 * 1.1, N = 9),
#    param2 = (name = "s", min = 0.958764 * 0.96, max = 0.958764 * 1.04, N = 9)
#    param2 = (name = "s", min = 0.1, max = 1.0, N = 9)
#    param2 = (name = "cosi", min = 0.28 * 0.9, max = 0.28 * 1.1, N = 9)
#    param2 = (name = "cosi", min = 0.0, max = 1.0, N = 9)
    param2 = (name = "gamma", min = 0.000773239 * 0.90, max = 0.000773239 * 1.1, N = 9)
#    param1 = (name = "m2", min = 0.0, max = 2.0, N = 9),
#    param2 = (name = "m1", min = 0.0, max = 2.0, N = 9)
    )

gsets = GridSetttings(
    N_refinement = 7,
    CL = [0.68],
    refinement_type = "massmass",
    delta_chisqr_max = 10.0,
    delta_chisqr_diff = 1.0,
    gr_in_chisqr = false
    )

pkf = PKFramework(test, obs_params, gsets)

calculate!(pkf, pf)
#-------------------------------------------------------------------------------------

grid_ddstg = load("/Users/abatrakov/Documents/PhD_work/projects/J2222-0137/saves/grid_J2222-0137_MPA1_nice_90CL.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=lvl_3σ), rasterized=true)
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[0.0, 0.5, 1.0, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1], colors="black")
#cs2 = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:m1], levels=[1.0, 1.5, 1.8], colors="black")
#clabel(cs2, cs2.levels)
#plot([], [], label=L"m_{\mathrm{p}}\,(M_\odot)", color="black")

cs = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[:chisqr] .- pkf.grid.params[:chisqr_min], levels=[lvl_68CL], colors="red")
#clabel(cs, cs.levels, fmt=Dict(lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\Delta\chi^{2}\, \, \mathrm{CL}", color="red")
ax.set_ylabel(get_label(pkf.test.param1.name), size=16)
ax.set_xlabel(get_label(pkf.test.param2.name), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)

cs3 = ax.contour(grid_ddstg.y, grid_ddstg.x, grid_ddstg.value[:chisqr] .- grid_ddstg.params[:chisqr_min], levels=levels=[lvl_90CL], colors=[cm(2)], linestyles=["-."])
plot([], [], label=L"\mathrm{Guo\ et\ al.}", color=cm(2), ls="-.", zorder=4)

#title(L"\mathrm{PK\ method, 23-26}")
legend(fontsize=12)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
#ax.invert_xaxis()
tight_layout()

#-------------------------------------------------------------------------------------

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
function get_label_PK(PK)
    if PK == :Pbdot
        return L"\dot{P}_\mathrm{b}"
    elseif PK == :k
        return L"\dot{\omega}"
    elseif PK == :gamma
        return L"\gamma"
    elseif PK == :varsigma
        return L"\varsigma"
    elseif PK == :h3
        return L"h_3"
    elseif PK == :r
        return L"r"
    elseif PK == :s
        return L"s"
    end
end

colors = [cm(0), cm(1), cm(2), cm(4), cm(8), cm(7)]
colors = [cm(9), cm(2), cm(4), cm(1), cm(8), cm(7)]

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
colors = [cm(1), cm(2), cm(6), cm(0), cm(3), cm(5)]

i = 1
for PK in keys(pkf.obs_params.PK)
    println((i,PK))
    if pkf.obs_params.PK[PK].err != 0.0
        PK_values = [pkf.obs_params.PK[PK].val - pkf.obs_params.PK[PK].err, pkf.obs_params.PK[PK].val, pkf.obs_params.PK[PK].val + pkf.obs_params.PK[PK].err]
        cs = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[PK], levels=PK_values, colors=[colors[i], colors[i], colors[i]], linestyles=[(0, (2, 2)), "-", (0, (2, 2))])
        plot([], [], label=get_label_PK(PK), color=colors[i], ls="-")
        i += 1
    end
end
legend(fontsize=12)
ax.set_xlabel(L"m_{\mathrm{p}},\,[M_\odot]", size=16)
ax.set_ylabel(L"m_{\mathrm{c}},\,[M_\odot]", size=16)
tight_layout()

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
function get_label_PK(PK)
    if PK == :Pbdot
        return L"\dot{P}_\mathrm{b}"
    elseif PK == :k
        return L"\dot{\omega}"
    elseif PK == :gamma
        return L"\gamma"
    elseif PK == :varsigma
        return L"\varsigma"
    elseif PK == :h3
        return L"h_3"
    elseif PK == :r
        return L"r"
    elseif PK == :s
        return L"s"
    end
end

function get_contour(grid, lvl, PK)
    cl = levels(contours(grid.x, grid.y, grid.value[PK], [lvl]))[1]
    X, Y = Float64[], Float64[]
    for line in lines(cl)
        Y_piece, X_piece = coordinates(line)
        Y = vcat(Y, Y_piece)
        X = vcat(X, X_piece)
    end
    XY = hcat(X, Y)
    XY = sortslices(XY,dims=1)
    return XY[:,1], XY[:,2]
end

colors = [cm(0), cm(1), cm(2), cm(4), cm(8), cm(7)]
colors = [cm(9), cm(2), cm(4), cm(1), cm(8), cm(7)]

cm = ColorMap(ColorSchemes.okabe_ito.colors, 8)
colors = [cm(1), cm(2), cm(6), cm(0), cm(3), cm(5)]
linestyles = ["-.", ":", "-", "--"]

i = 1
for PK in keys(pkf.obs_params.PK)
    println((i,PK))
    if pkf.obs_params.PK[PK].err != 0.0
        PK_values = [pkf.obs_params.PK[PK].val - pkf.obs_params.PK[PK].err, pkf.obs_params.PK[PK].val + pkf.obs_params.PK[PK].err]
        cs = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[PK], levels=PK_values, colors=[colors[i], colors[i]], linestyles=[linestyles[i], linestyles[i]], linewidth=1.0)
        plot([], [], label=get_label_PK(PK), color=colors[i], ls=linestyles[i])
        itp1 = LinearInterpolation(get_contour(pkf.grid, PK_values[1], PK)..., extrapolation_bc=Line())
        itp2 = LinearInterpolation(get_contour(pkf.grid, PK_values[2], PK)..., extrapolation_bc=Line())
        fill_between(pkf.grid.y, itp1.(pkf.grid.y), itp2.(pkf.grid.y), color=colors[i], alpha=0.2)
        if PK == :Pbdot
            cs = ax.contour(pkf.grid.y, pkf.grid.x, pkf.grid.value[PK], levels=[pkf.obs_params.PK[PK].val], colors=[colors[i]], linestyles=["-"])
        end
        i += 1
    end
end
ylim(pkf.grid.x[1], pkf.grid.x[end])
legend(fontsize=12)
ax.set_xlabel(L"m_{\mathrm{p}},\,[M_\odot]", size=16)
ax.set_ylabel(L"m_{\mathrm{c}},\,[M_\odot]", size=16)
tight_layout()

#-------------------------------------------------------------------------------------
