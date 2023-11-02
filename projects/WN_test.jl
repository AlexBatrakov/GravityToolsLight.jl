using Revise
using PyPlot
using Contour
using JLD
using ColorSchemes
using Statistics
using DelimitedFiles
using Distributions
using NLsolve
using Optim
using Roots


cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN")
cd("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/model_WN")
#tim_file = readdlm("J1141-6545_until_2018.tim", String)
#tim_file = readdlm("UWL_partial.tim", String)
tim_file = readdlm("J1141-6545_full.tim", String)
N_TOAs = size(tim_file)[1] - 2

postfit = readdlm("postfit.res", Float64)[:,2]
unc = [parse(Float64, tim_file[i,4]) for i in 3:length(tim_file[:,3])]
trans(unc, EFAC, EQUAD) = sqrt(EFAC^2*unc^2 + EQUAD^2)
res = postfit .* 1e6

sqrt(mean(res .^ 2))
mean((res ./ unc) .^2)

unc_tr = copy(unc)


mean_fun(res, unc, EFAC, EQUAD) = (EFAC == EQUAD == 0) ? Inf : mean(res ./ trans.(unc, EFAC, EQUAD) ) # -> 0.0
std_fun(res, unc, EFAC, EQUAD) = (EFAC == EQUAD == 0) ? Inf : std(res ./ trans.(unc, EFAC, EQUAD) )   # -> 1.0
skew_fun(res, unc, EFAC, EQUAD) = (EFAC == EQUAD == 0) ? Inf : skewness(res ./ trans.(unc, EFAC, EQUAD) ) # -> 0.0
kurt_fun(res, unc, EFAC, EQUAD) = (EFAC == EQUAD == 0) ? Inf : kurtosis(res ./ trans.(unc, EFAC, EQUAD) )  # -> 0.0 
#chisqr_fun(res, unc, EFAC, EQUAD) = mean(res .^ 2 ./ trans.(unc, EFAC, EQUAD) .^2 )

EFAC_arr = collect(0:0.01:5)
EQUAD_arr = collect(0:1:300)

#std_arr = [std_fun(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]
#kurt_arr = [kurt_fun(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

#minimum(kurt_arr)

# fig, ax = subplots()
# pclm = ax.pcolormesh(EQUAD_arr, EFAC_arr, kurt_arr, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=3.0))
# cbar = colorbar(pclm)
# cs_std = ax.contour(EQUAD_arr, EFAC_arr, std_arr, levels = [1.0], colors=["green"])
# cs_kurt = ax.contour(EQUAD_arr, EFAC_arr, kurt_arr, levels = [0.0], colors=["red"])


be_arr = ["CASPSR", "PDFB1", "PDFB2", "PDFB3", "MOPSR", "SCAMP", "CPSR2_LOW", "CPSR2_HIGH", "Medusa"]

cut_sigma_level = 4.0

inds_dict = Dict()
inds_cut_dict = Dict()
inds_length_dict = Dict()
mean_dict = Dict()
std_dict = Dict()
skew_dict = Dict()
kurt_dict = Dict()
EFAC_dict = Dict()
EQUAD_dict = Dict()
log10EQUAD_dict = Dict()

function estimate_WN(res, unc, inds, inds_cut)
    lambda_fun(lambda) = kurt_fun(res[inds_cut], unc[inds_cut], 1.0, lambda)
    sol_lambda = optimize(lambda_fun, [0.0])
    lambda = Optim.minimizer(sol_lambda)[1]
    lambda = lambda > 0.0 ? lambda : 0.0
    EFAC_fun(EFAC) = std_fun(res[inds], unc[inds], EFAC, EFAC*lambda) - 1.0
    EFAC = find_zero(EFAC_fun, 1.0)
    EQUAD = lambda * EFAC
    return EFAC, EQUAD
end

estimate_WN(be) = estimate_WN(res, unc, inds_dict[be], inds_cut_dict[be])

for be in be_arr
    inds_dict[be] = inds = filter(ind -> be in tim_file[ind+2, :], collect(1:N_TOAs))
    inds_cut = copy(inds)
    
    for i in 1:10
        cut_level = std(res[inds_cut] ./ unc[inds_cut]) * cut_sigma_level
        inds_cut = filter(ind -> abs(res[ind] / unc[ind]) < cut_level, inds_cut)
    end
    inds_cut_dict[be] = inds_cut

    EFAC_dict[be], EQUAD_dict[be] = estimate_WN(be)
    log10EQUAD_dict[be] = log10(EQUAD_dict[be] > 0 ? EQUAD_dict[be] : NaN) - 6.0

    unc_tr[inds] = trans.(unc[inds], EFAC_dict[be], EQUAD_dict[be])

    inds_length_dict[be] = length(inds), length(inds_cut)
    mean_dict[be]        = mean(res[inds] ./ unc[inds]),     mean(res[inds_cut] ./ unc[inds_cut]),     mean(res[inds] ./ unc_tr[inds]),     mean(res[inds_cut] ./ unc_tr[inds_cut])
    std_dict[be]         = std(res[inds] ./ unc[inds]),      std(res[inds_cut] ./ unc[inds_cut]),      std(res[inds] ./ unc_tr[inds]),      std(res[inds_cut] ./ unc_tr[inds_cut])
    skew_dict[be]        = skewness(res[inds] ./ unc[inds]), skewness(res[inds_cut] ./ unc[inds_cut]), skewness(res[inds] ./ unc_tr[inds]), skewness(res[inds_cut] ./ unc_tr[inds_cut])
    kurt_dict[be]        = kurtosis(res[inds] ./ unc[inds]), kurtosis(res[inds_cut] ./ unc[inds_cut]), kurtosis(res[inds] ./ unc_tr[inds]), kurtosis(res[inds_cut] ./ unc_tr[inds_cut])
end

inds_length_dict
mean_dict
std_dict
skew_dict
kurt_dict
EFAC_dict
EQUAD_dict
log10EQUAD_dict

function print_WN_params(be_arr, EFAC_dict, log10EQUAD_dict)
    for be in be_arr
            println("TNEF -be $be ", EFAC_dict[be])
    end
    for be in be_arr
        if !isnan(log10EQUAD_dict[be])
            println("TNEQ -be $be ", log10EQUAD_dict[be])
        end
    end
end

print_WN_params(be_arr, EFAC_dict, log10EQUAD_dict)

#-------------------------------------------------------------------------------------

be = "MOPSR"
inds = inds_dict[be]
inds_cut = inds_cut_dict[be]


std_arr = [std_fun(res[inds], unc[inds], EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]
kurt_arr = [kurt_fun(res[inds], unc[inds], EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

minimum(kurt_arr)

fig, ax = subplots()
pclm = ax.pcolormesh(EQUAD_arr, EFAC_arr, kurt_arr, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=10.0))
cbar = colorbar(pclm)
cs_std = ax.contour(EQUAD_arr, EFAC_arr, std_arr, levels = [1.0], colors=["green"])
cs_kurt = ax.contour(EQUAD_arr, EFAC_arr, kurt_arr, levels = sort([minimum(kurt_arr)*1.0001, minimum(kurt_arr)*0.9999]), colors=["red","pink"])
#cs_mean = ax.contour(EQUAD_arr, EFAC_arr, mean_arr, levels = [0.0], colors=["black"])


#-------------------------------------------------------------------------------------


EFAC_dict = Dict("CASPSR" => 1.0358090018551152,
       "PDFB1" => 1.2424153969446858,
       "PDFB2" => 1.0750588236136025,
       "PDFB3" => 1.3267880918752266,
       "MOPSR" => 0.9871760803080492,
       "SCAMP" => 1.3617029106000096,
       "CPSR2_LOW" => 1.2086765660697325,
       "CPSR2_HIGH" => 1.2050174255728376,
       "Medusa" => 1.2083866324857806)

EQUAD_dict = Dict("CASPSR" => exp10(6 + -4.249931277939291),
       "PDFB1" => exp10(6 + -4.268684581358633),
       "PDFB2" => exp10(6 + -4.198238266972263),
       "PDFB3" => exp10(6 + -4.400461882182785),
       "MOPSR" => exp10(6 + -4.55704095459974),
       "SCAMP" => exp10(6 + -3.990612961320208),
       "CPSR2_LOW" => exp10(6 + -4.291451748974723),
       "CPSR2_HIGH" => exp10(6 + -4.298204380047417),
       "Medusa" => exp10(6 + -4.378062085685085))