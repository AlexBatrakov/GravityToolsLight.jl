using Revise
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
using DelimitedFiles
using Distributions
pygui(true)

using Distributed
addprocs(8)
#pkill -f /opt/julia-1.8.2/bin/julia

@everywhere using GravityToolsLight

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

#-------------------------------------------------------------------------------------

cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN")
# cd /Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_and_12k_gauss
# tempo2 -f J1141-6545_until_2018_DDSTG.par J1141-6545_until_2018.tim -nobs 21000 -newpar J1141-6545_until_2018_DDSTG_out.par -writeres

prefit = readdlm("prefit.res")
postfit = readdlm("postfit.res")
residuals = readdlm("residuals.dat")

#plot(prefit[:,1], prefit[:,2], ".")
#errorbar(prefit[:,1], prefit[:,2], yerr = prefit[:,3], fmt =".")

tim_file = readdlm("J1141-6545_until_2018.tim", String)
tim_file = readdlm("J1141-6545_until_2018_pure_GR.tim", String)
#tim_file = readdlm("J1141-6545_full_shift_pure_GR.tim", String)
#tim_file = readdlm("J1141-6545_full_pure_GR4.tim", String)
#tim_file = readdlm("J1141-6545_until_2018_pure_GR3.tim", String)
tim_file = readdlm("UWL_partial.tim", String)
tim_file = readdlm("UWL_partial_pure_GR.tim", String)
tim_file = readdlm("J1141-6545_full.tim", String)
tim_file = readdlm("J1141-6545_full_pure_GR.tim", String)
#tim_file = readdlm("J1141-6545_until_2018_WN2.tim", String)

TOA_WN = randn(length(tim_file[:,3]))

for i in 3:length(tim_file[:,3])
    TOA_int_string, TOA_float_string = split(tim_file[i,3], ".")
#    TOA_err = parse(Float64, tim_file[i,4])
    TOA_int = parse(Int64, TOA_int_string)
    TOA_float = parse(Float64, "." * TOA_float_string)

#    TOA_res = prefit[i-2,2] .- postfit[i-2,2]
    TOA_res = residuals[i-2,2]

    TOA_float = TOA_float - 0.0 * TOA_res / (24 * 3600) + 1.0 * TOA_err[i-2] * TOA_WN[i-2] * 1e-6 / (24 * 3600)
    if TOA_float > 1
        TOA_int += 1
        TOA_float -= 1.0
    elseif TOA_float < 0
        TOA_int -= 1
        TOA_float += 1.0
    end

    TOA_float += 1.0

    tim_file[i,3] = string(TOA_int) * string(TOA_float)[2:end]
end


open("J1141-6545_full_WN.tim", "w") do io
#open("J1141-6545_full_shift_antiWN.tim", "w") do io
#open("J1141-6545_until_2018_WN.tim", "w") do io
#open("UWL_partial_pure_GR.tim", "w") do io
#open("testUWL.tim", "w") do io
#open("J1141-6545_until_2018_pure_GR.tim", "w") do io
    writedlm(io, tim_file, ' ')
end

prefit = readdlm("prefit.res")
postfit = readdlm("postfit.res")
residuals = readdlm("residuals.dat")

fig, ax = subplots()
plot(residuals[:,1] .+ 54000.0, postfit[:,2], ".", label="postfit")
plot(residuals[:,1] .+ 54000.0, residuals[:,2], ".", label="resid")
plot(residuals[:,1] .+ 54000.0, prefit[:,2] .- postfit[:,2], ".", label="prefit - postfit")
legend()
title("J1141-6545_full_pure_GR5.tim")
savefig("residuals_pure_GR5.pdf")


res = residuals[:,2] .* 1e6
unc = [parse(Float64, tim_file[i,4]) for i in 3:length(tim_file[:,3])]

fig, ax = subplots()
scatter(unc, abs.(res), c=res.^2 ./ (unc) .^2, norm = matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e+3), cmap="jet")
colorbar()
yscale("log")
xscale("log")
xlim(left=5e0, right=2e3)
ylim(bottom=2e-3, top=3.5e3)


u_arr = 10.0 .^ collect(LinRange(0,4,101))
plot(u_arr, 1.0 * u_arr .+ 0.0, "black")

out = readdlm("test.out")



tim_file = readdlm("J1141-6545_until_2018.tim", String)
prefit = readdlm("prefit.res")
postfit = readdlm("postfit.res")
residuals = readdlm("residuals.dat")

#res = residuals[:,2] .* 1e6
res = postfit .* 1e6
unc = [parse(Float64, tim_file[i,4]) for i in 3:length(tim_file[:,3])]

A_arr = collect(0:0.01:3)
B_arr = collect(0:1:300)

chisqr_fun(A,B, res=res, unc=unc) = mean( (res .^ 2 ./ ( (A*unc) .^2  .+ B .^ 2  )) )
test_fun(A,B, res=res, unc=unc) = sqrt(mean( (abs.(res) .- sqrt.( (A*unc) .^2  .+ B .^2)) .^2 ))

test_arr = [test_fun(A, B) for A in A_arr, B in B_arr]
chisqr_arr = [chisqr_fun(A, B) for A in A_arr, B in B_arr]

fig, ax = subplots()
pclm = ax.pcolormesh(B_arr, A_arr, test_arr .- minimum(test_arr), cmap="Blues_r")
cbar = colorbar(pclm)

cs = ax.contour(B_arr, A_arr, chisqr_arr, levels = [1.0], colors=["red"])
cs = ax.contour(B_arr, A_arr, test_arr .- minimum(test_arr), levels = [0.5, 1.0, 2.0, 3.0, 4.0, 50.0], colors=["green"])

#-------------------------------------------------------------------------------------
using DelimitedFiles
using Statistics
using PyPlot
using Distributions

cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN")
tim_file = readdlm("J1141-6545_until_2018.tim", String)
#tim_file = readdlm("UWL_partial.tim", String)
postfit = readdlm("postfit.res", Float64)[:,2]
unc = [parse(Float64, tim_file[i,4]) for i in 3:length(tim_file[:,3])]
trans(unc, A=1.5, B=30.0 ) = sqrt(A^2*unc^2 + B^2)
#res_tr = [trans.(unc)*randn() for unc in unc]
res = postfit .* 1e6

std_fun(A,B) = std(res ./ trans.(unc, A,B) )   # -> 1.0
kurt_fun(A,B) = kurtosis(res ./ trans.(unc, A,B) )  # -> 0.0
#fun2(A,B) = mean( (res_tr ./ trans.(unc, A,B)) .^2 ) # -> 1.0

A_arr = collect(0:0.01:3)
B_arr = collect(0:1:300)

std_arr = [std_fun(A, B) for A in A_arr, B in B_arr]
kurt_arr = [kurt_fun(A, B) for A in A_arr, B in B_arr]

kurt_arr[1,1] = Inf
minimum(kurt_arr)

fig, ax = subplots()
pclm = ax.pcolormesh(B_arr, A_arr, kurt_arr, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=3.0))
cbar = colorbar(pclm)
cs_std = ax.contour(B_arr, A_arr, std_arr, levels = [1.0], colors=["green"])
cs_kurt = ax.contour(B_arr, A_arr, kurt_arr, levels = [0.018], colors=["red"])


be_dict = ["CASPSR", "PDFB1", "PDFB2", "PDFB3", "MOPSR", "SCAMP", "CPSR2_LOW", "CPSR2_HIGH"]

inds_dict = Dict()

for be in be_dict
    inds_dict[be] = inds = filter(ind -> be in tim_file[ind, :], collect(1:20861))
end

be = "SCAMP"
inds = inds_dict[be]

mean(res[inds] ./ unc[inds])
std(res[inds] ./ unc[inds])
skewness(res[inds] ./ unc[inds])
kurtosis(res[inds] ./ unc[inds])

for i in 1:10
    cut = std(res[inds] ./ unc[inds]) * 3
    inds = filter(ind -> abs(res[ind] / unc[ind]) < cut, inds)
    std(res[inds] ./ unc[inds])
end

mean_fun(A,B, inds=inds) = mean(res[inds] ./ trans.(unc[inds], A,B) )
std_fun(A,B, inds=inds) = std(res[inds] ./ trans.(unc[inds], A,B) )
skew_fun(A,B, inds=inds) = skewness(res[inds] ./ trans.(unc[inds], A,B) )
kurt_fun(A,B, inds=inds) = kurtosis(res[inds] ./ trans.(unc[inds], A,B) ) 
chisqr_fun(A,B, inds=inds) = mean((res[inds] .^ 2 ./ trans.(unc[inds], A,B) .^2) )
 
 
A_arr = collect(0:0.002:5)
B_arr = collect(0:0.2:300)

#mean_arr = [mean_fun(A, B) for A in A_arr, B in B_arr]
std_arr = [std_fun(A, B) for A in A_arr, B in B_arr]
#skew_arr = [skew_fun(A, B) for A in A_arr, B in B_arr]
kurt_arr = [kurt_fun(A, B) for A in A_arr, B in B_arr]
#chisqr_arr = [chisqr_fun(A, B) for A in A_arr, B in B_arr]

#skew_arr[1,1] = Inf
#mean_arr[1,1] = Inf
#minimum(mean_arr)
kurt_arr[1,1] = Inf
minimum(kurt_arr)

fig, ax = subplots()
pclm = ax.pcolormesh(B_arr, A_arr, kurt_arr, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0))
cbar = colorbar(pclm)
cs_std = ax.contour(B_arr, A_arr, std_arr, levels = [1.0], colors=["green"])
cs_kurt = ax.contour(B_arr, A_arr, kurt_arr, levels = sort([minimum(kurt_arr)*1.0001, minimum(kurt_arr)*0.9999]), colors=["red","pink"])
#cs_mean = ax.contour(B_arr, A_arr, mean_arr, levels = [0.0], colors=["black"])


#-------------------------------------------------------------------------------------

cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN")

#tim_file = readdlm("J1141-6545_until_2018.tim", String)
tim_file = readdlm("UWL_partial.tim", String)
unc_RN = [parse(Float64, tim_file[i,4]) for i in 3:length(tim_file[:,3])]

#out_file = readdlm("until_2018.out", Float64)
out_file = readdlm("UWL_partial_F5.out", Float64)
resid = readdlm("residuals.dat", Float64)[:,2]
postfit = readdlm("postfit.res", Float64)[:,2]

freq = out_file[:,2]
unc_TN = out_file[:,3]
prefit = out_file[:,4] .* 1e6
#postfit = out_file[:,5] .* 1e6

res_TN = postfit .* 1e6

fig, ax = subplots()
scatter(unc_RN, abs.(res_TN), c=res_TN.^2 ./ (unc_RN) .^2, norm = matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e+3), cmap="jet")
colorbar()
yscale("log")
xscale("log")
xlim(left=5e0, right=2e3)
ylim(bottom=2e-3, top=5e3)

fig, ax = subplots()
scatter(unc_TN, abs.(res_TN), c=res_TN.^2 ./ (unc_TN) .^2, norm = matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e+3), cmap="jet")
colorbar()
yscale("log")
xscale("log")
xlim(left=5e0, right=2e3)
ylim(bottom=2e-3, top=5e3)

std(res_TN ./ unc_RN)
std(res_TN ./ unc_TN)

kurtosis(res_TN ./ unc_RN)
kurtosis(res_TN ./ unc_TN)

be_dict = ["CASPSR", "PDFB1", "PDFB2", "PDFB3", "MOPSR", "SCAMP", "CPSR2_LOW", "CPSR2_HIGH"]

inds_dict = Dict()
std_dict = Dict()
kurt_dict = Dict()


for be in be_dict
    inds_dict[be] = inds = filter(ind -> be in tim_file[ind, :], collect(1:20861))
    std_dict[be] = std(res_TN[inds] ./ unc_TN[inds])
    kurt_dict[be] = kurtosis(res_TN[inds] ./ unc_TN[inds])
end

n = 0
for be in be_dict
    n += length(inds_dict[be])
end
n

inds = inds_dict["PDFB1"]
unc = unc_TN
res = res_TN
fig, ax = subplots()
scatter(unc[inds], abs.(res[inds]), c=res[inds].^2 ./ (unc[inds]) .^2, norm = matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e+3), cmap="jet")
colorbar()
yscale("log")
xscale("log")
xlim(left=5e0, right=2e3)
ylim(bottom=2e-3, top=5e3)

fig, ax = subplots()
plot(unc_RN[inds], unc_TN[inds], ".")

fig, ax = subplots()
scatter(unc_RN[inds], unc_TN[inds], c = freq[inds])
colorbar()
yscale("log")
xscale("log")

#-------------------------------------------------------------------------------------
# cd /software//tempo2/T2runtime
# ln -s /Users/abatrakov/Documents/Work/PhD/computed_grids_fine data_ddstg
# cd /Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/
# cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/fake_with_gauss")
cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN")

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

# tsets = TempoSettings(
#     par_file_init = "J1141-6545_until_2018_DDSTG.par",
#     tim_file = "J1141-6545_until_2018.tim",
#     add_flag = "-nobs 21000 -newpar",
#     fit_XPBDOT = false,
#     nits_first_step = 3,
#     gain_fisrt_step = 1.0
#     )

tsets = TempoSettings(
    par_file_init = "J1141-6545_full_DDSTG_best2_WN.par",
#    tim_file = "J1141-6545_until_2018_WN.tim",
#    tim_file = "J1141-6545_full_WN.tim",
#    tim_file = "J1141-6545_until_2018_shift_antiWN.tim",
#    tim_file = "J1141-6545_full_shift_antiWN.tim",
    tim_file = "J1141-6545_until_2018_WN.tim",
#    tim_file = "UWL_partial_WN.tim",
#    add_flag = "-nobs 34000 -newpar -set XPBDOT -3.4536e-15",
    add_flag = "-nobs 34000 -newpar",
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


rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()
pclm = ax.pcolormesh(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr] .- tf.grid.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=tf.gsets.delta_chisqr_max), rasterized=true)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")
cs = ax.contour(tf.grid.y, tf.grid.x, tf.grid.value[:chisqr_cut] .- tf.grid.params[:chisqr_min], levels=tf.gsets.contours, colors=["red"])

save("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_full_WN_BSk22.jld",  "grid", tf.grid)


#-------------------------------------------------------------------------------------

grid_full = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_full_WN_BSk22.jld", "grid")
grid_full_x = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_full_WN_XPBDOT_BSk22.jld", "grid")
grid_ntil_2018 = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_until_2018_WN_BSk22.jld", "grid")

#-------------------------------------------------------------------------------------

grid_DDSTG_J1141_6545_full_shift_antiWN = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_full_shift_antiWN_BSk22.jld", "grid")
grid_DDSTG_J1141_6545_until_2018_shift_antiWN = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_DDSTG_J1141-6545_until_2018_shift_antiWN_BSk22.jld", "grid")
grid_PK_Triple_System = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Triple_System_BSk22.jld", "grid")
grid_PK_Double_Pulsar = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_Double_Pulsar_BSk22.jld", "grid")
grid_PK_J1738_0333 =   load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/saves/grid_PK_J1738+0333_Guo_BSk22.jld", "grid")

rc("mathtext",fontset="cm")
rc("font", family="serif", size=12)
cm = ColorMap(ColorSchemes.seaborn_colorblind.colors, 10)
fig, ax = subplots()

grid_bg = grid_full 

pclm = ax.pcolormesh(grid_bg.y, grid_bg.x, grid_bg.value[:chisqr] .- grid_bg.params[:chisqr_min], cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=6.0), rasterized=true)
cbar = colorbar(pclm)
cbar.set_label(L"$\Delta\chi^{2}$")

cs_bg = ax.contour(grid_bg.y, grid_bg.x, grid_bg.value[:chisqr_cut] .- grid_bg.params[:chisqr_min], levels=[lvl_90CL], colors=["red"])
plot([], [], label=L"\mathrm{DDSTG, WN, 22 yr}", color="red", ls="-")
#plot_contour(ax, grid_DDSTG_J1141_6545_full_shift_antiWN,  lvl_90CL, color=cm(3), label=L"\mathrm{J}1141-6545\ \mathrm{DDSTG, full}", N_smooth=10)

grid2 = grid_full_x

cs2 = ax.contour(grid2.y, grid2.x, grid2.value[:chisqr_cut] .- grid2.params[:chisqr_min], levels=[lvl_90CL], colors=["red"], linestyles=["--"])
plot([], [], label=L"\mathrm{DDSTG, WN, 22 ye, XPBDOT}", color="red", ls="--")

grid3 = grid_until_2018

cs2 = ax.contour(grid3.y, grid3.x, grid3.value[:chisqr_cut] .- grid3.params[:chisqr_min], levels=[lvl_90CL], colors=["blue"], linestyles=["-"])
plot([], [], label=L"\mathrm{DDSTG, WN, 18 yr}", color="blue", ls="--")


grid_DDSTG_J1141_6545 = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/grid_DDSTG_J1141_6545_BSk22.jld", "grid")

cs_DDSTG_J1141_6545 = ax.contour(grid_DDSTG_J1141_6545.y, grid_DDSTG_J1141_6545.x, grid_DDSTG_J1141_6545.value[:chisqr] .- grid_DDSTG_J1141_6545.params[:chisqr_min], levels=[lvl_90CL], colors=["red"], linestyles=[":"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{DDSTG, RN, 18 yr}", color="red", ls=":")



cs_Triple_System = ax.contour(grid_PK_Triple_System.y, grid_PK_Triple_System.x, grid_PK_Triple_System.value[:chisqr] .- grid_PK_Triple_System.params[:chisqr_min], levels=[lvl_90CL], colors=["green"], linestyles=["-"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{Triple System}", color="green", ls="-")

cs_Double_Pulsar = ax.contour(grid_PK_Double_Pulsar.y, grid_PK_Double_Pulsar.x, grid_PK_Double_Pulsar.value[:chisqr] .- grid_PK_Double_Pulsar.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(4)], linestyles=["-."])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}0737-3039\mathrm{A}}", color=cm(4), ls="-.")

cs_J1738_0333 = ax.contour(grid_PK_J1738_0333.y, grid_PK_J1738_0333.x, grid_PK_J1738_0333.value[:chisqr] .- grid_PK_J1738_0333.params[:chisqr_min], levels=[lvl_90CL], colors=[cm(2)], linestyles=[":"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{J}1738+0333", color=cm(2), ls=":")


ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
title("BSk22")
#ax.invert_xaxis()
xlim(left=-6.0, right=+10.0)
ylim(top=-2.0, bottom=-4.0)
tight_layout()

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

plot_contour(ax, grid_DDSTG_J1141_6545,  lvl_90CL, color=cm(3), label=L"Science data", N_smooth=10)

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

cs_Triple_System = ax.contour(grid_PK_Triple_System.y, grid_PK_Triple_System.x, grid_PK_Triple_System.value[:chisqr] .- grid_PK_Triple_System.params[:chisqr_min], levels=[lvl_90CL], colors=["red"], linestyles=["-"])
#clabel(cs, cs.levels, fmt=Dict(1.0 => "1.0", lvl_68CL => "σ", lvl_95CL => "2σ"))
plot([], [], label=L"\mathrm{Triple System}", color="red", ls="-")


ax.set_ylabel(get_label("log10alpha0"), size=16)
ax.set_xlabel(get_label("beta0"), size=16)
ax.axhline(y=log10(3.4e-3), color=cm(7), label=L"\mathrm{Cassini}", ls="--", zorder=5)
legend(fontsize=12)
#ax.invert_xaxis()
xlim(left=-6.0, right=+10.0)
ylim(top=-1.5, bottom=-4.0)
tight_layout()

savefig("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/science_paper/plot_BSk22.pdf")