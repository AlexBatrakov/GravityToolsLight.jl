using Revise
using GravityToolsLight
using PyPlot
using HypothesisTests
using DelimitedFiles
using Distributions
pygui(true)


basic_settings = bsets = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3), TP("TNRedGam", 2.5), TP("TNRedC", 33)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=true)
    )


# results_basic = run_tempo_basic(basic_settings)

global_iters_settings = gisets = GlobalIterationsSettings(
    keys = GlobalIterationsKeys(iterative_mode=true, save_global_iterations=true),
    iters = 4,
    nits = [2,1,1,1],
    gain = [1.0, 1.0, 1.0, 1.0]
    )

# results_global_iters = run_tempo_global_iters(basic_settings, global_iters_settings)

parameter_sweep_settings = ParameterSweepSettings(
    parameter_name = "TNRedAmp",
    values = collect(-12:0.5:-7)
)

sweep_results = run_tempo_parameter_sweep(basic_settings, global_iters_settings, parameter_sweep_settings)

chisqr_arr_2p5 = [sweep_results[iter].last_internal_iteration.result.chisqr for iter in 1:length(parameter_sweep_settings.values)]
TRES_arr_2p5 = [sweep_results[iter].last_internal_iteration.result.TRES.post_fit for iter in 1:length(parameter_sweep_settings.values)]
plot(parameter_sweep_settings, TRES_arr)

chisqr_arr
TRES_arr

settings = basic_settings

tim_file_path = joinpath(settings.work_dir, settings.tim_file)
backends = settings.backends

# Чтение данных из файла .tim
tim_file_data = readdlm(tim_file_path, String)
uncertainties = parse.(Float64, tim_file_data[3:end, 4])
residuals = readdlm(joinpath(settings.work_dir, "postfit.res"), Float64)[:, 2] .* 1e6


EFAC_max = 3.0
EQUAD_max = 200.0
lambda_max = 1000.0

EFAC_arr = LinRange(0, EFAC_max, 501)
EQUAD_arr = LinRange(0, EQUAD_max, 501)
lambda_arr = LinRange(0, lambda_max, 501)



#-------------------------------------------------------------------------------------




EFACs, EQUADs, log10EQUADs = calculate_EFACs_EQUADs(basic_settings)

transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

res = residuals[indices]
unc = uncertainties[indices]

fig, ax = subplots()
hist(res, bins = 100)

fig, ax = subplots()
hist(unc, bins = 100, alpha=0.5)
hist(abs.(res), bins = 100, alpha=0.5)

fig, ax = subplots()
hist(abs.(res ./ unc), bins = 100, alpha=0.5)
hist(abs.(randn(length(res))), bins = 100, alpha=0.5)

N = length(res)
EFAC, lambda = 1.0, 84
EQUAD = EFAC * lambda

unc_tr = transform_uncertainty.(unc, EFAC, EQUAD)

mean(abs.(res ./ unc)), mean(abs.(res ./ unc_tr)), sqrt(2.0 / pi), mean(abs.(randn(N)))
std(abs.(res ./ unc)), std(abs.(res ./ unc_tr)), sqrt(1.0 - 2.0 / pi), std(abs.(randn(N)))
skewness(abs.(res ./ unc)), skewness(abs.(res ./ unc_tr)), sqrt(2.0) * (4.0 - pi) / (pi - 2.0) ^ 1.5, skewness(abs.(randn(N)))
kurtosis(abs.(res ./ unc)), kurtosis(abs.(res ./ unc_tr)), 8.0 * (pi - 3.0) / (pi - 2.0) ^ 2 , kurtosis(abs.(randn(N)))

fig, ax = subplots()
counts, bins = hist(abs.(res ./ unc_tr), bins = 100, alpha=0.5, density = true)
counts_n, bins = hist(abs.(randn(length(res) * 100)), bins = bins, alpha=0.5, density = true)

fig, ax = subplots()
plot(bins[1:end-1], counts ./ counts_n)
axhline(y=1.0, color="red")

Q = filter(q -> !isnan(q) && !isinf(q) , counts ./ counts_n)
mean(Q)
std(Q)

res_std = res ./ std(res)
unc_std = unc ./ std(unc)

trans(unc, EFAC, EQUAD) = transform_uncertainty.(unc, EFAC, EQUAD) ./ std(transform_uncertainty.(unc, EFAC, EQUAD))


plot(res_std, res_std ./ unc_std, ".")
plot(res_std, res_std ./ trans(unc, 0.0, 100.0), ".")
plot(res_std, res_std ./ trans(unc, 4.0, 0.0), ".")


alpha(res_std, unc, EFAC, EQUAD) = atan.(res_std ./ trans(unc, EFAC, EQUAD) ./     std(residuals[indices] ./ uncertainties[indices]) ./ uncertainties[indices])

plot(atan.(std(residuals[indices] ./ uncertainties[indices]) ./ uncertainties[indices]), residuals[indices] ./ uncertainties[indices] ./ std(residuals[indices] ./ uncertainties[indices]), ".")


plot(unc_std, res_std, ".")

plot(res_std ./ unc_std, atan.(res_std ./ unc_std, res_std),".")


#-------------------------------------------------------------------------------------



EFAC_max = 3.0
EQUAD_max = 150.0
lambda_max = 500.0

EFAC_arr = LinRange(0, EFAC_max, 1001)
EQUAD_arr = LinRange(0, EQUAD_max, 1001)
lambda_arr = LinRange(0, lambda_max, 1001)

KS_objective(res, unc, EFAC, EQUAD/EFAC)

adaptive_ks_fun = AdaptiveSparseGrid([0.0, 0.0], [EFAC_max, lambda_max], max_depth = 10, tol=1e-4) do (EFAC, lambda)
    res_norm = res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)
    ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
    return (ks_stat = KS_objective(res, unc, EFAC, lambda), 
    )
    end


ks_arr = [EFAC == 0 ? Inf : adaptive_ks_fun([EFAC, lambda]) for EFAC in EFAC_arr, lambda in lambda_arr]

EFAC, EQUAD = GravityToolsLight.estimate_WhiteNoise_KS(residuals, uncertainties, indices)
lambda = EQUAD / EFAC

delete_nans(x) = isnan(x) ? Inf : x

ks_arr .= delete_nans.(ks_arr)

fig, ax = subplots()
pcolormesh(lambda_arr, EFAC_arr, log10.(ks_arr))
colorbar()
plot(lambda, EFAC, "x", color="red")


#ks_arr_true = [EFAC == 0 ? Inf : KS_objective(res, unc, EFAC, lambda) for EFAC in EFAC_arr, lambda in lambda_arr]
function get_best_ks_lambda(lambda)
    try
    EFAC_function(EFAC) = std(res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)) - 1.0
    EFAC_init = find_zero(EFAC_function, 1.0)
    #print(lambda, " ", EFAC_init)

    ks_fun_for_EFAC(EFAC) = KS_objective(res, unc, EFAC, lambda)
    EFAC_opt = optimize(ks_fun_for_EFAC, [EFAC_init])
    EFAC_best = Optim.minimizer(EFAC_opt)[1]
    #println(" ", EFAC_best)
    return ks_fun_for_EFAC(EFAC_best)
    catch error
        return Inf
    end
end


ks_min_lambda = [minimum(ks_arr[:,i]) for i in eachindex(lambda_arr)]
ks_min_lambda_true = [get_best_ks_lambda(lambda) for lambda in lambda_arr]
#ks_min_lambda_true = [minimum(ks_arr_true[:,i]) for i in eachindex(lambda_arr)]
fig, ax = subplots()
plot(lambda_arr, ks_min_lambda, label = "interp")
plot(lambda_arr, ks_min_lambda_true, label = "true")
axvline(x=lambda, color="red")
legend()



function get_best_ks_EFAC(EFAC)
    try
        lambda_function(lambda) = std(res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)) - 1.0
        lambda_init = find_zero(lambda_function, 1)

    ks_fun_for_EFAC(lambda) = KS_objective(res, unc, EFAC, lambda)
    lambda_opt = optimize(ks_fun_for_EFAC, [lambda_init])
    lambda_best = Optim.minimizer(lambda_opt)[1]
    #println(" ", EFAC_best)
    return ks_fun_for_EFAC(lambda_best)
    catch error
        return Inf
    end
end

ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]
ks_min_EFAC_true = [get_best_ks_EFAC(EFAC) for EFAC in EFAC_arr]
#ks_min_EFAC_true = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]
fig, ax = subplots()
plot(EFAC_arr, ks_min_EFAC, label = "interp")
plot(EFAC_arr, ks_min_EFAC_true, label = "true")
axvline(x=EFAC, color="red")
legend()



#-------------------------------------------------------------------------------------
basic_settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/full_WN",
    version = Tempo2(),
    par_file_init = "J1141-6545_until_2018_DDSTG.par",
    tim_file = "J1141-6545_until_2018.tim",
    flags = "-nobs 34000 -newpar -writeres -residuals",
    tparams = [TP("EOS", "BSk22"), TP("COMP_TYPE", "WD"), TP("ALPHA0", 0.0), TP("BETA0", 0.0), TP("NITS", 3)],
    keys = BasicTempoKeys(silent=true, print_output=true, save_internal_iterations=true, fit_EFACs_EQUADs=false)
    )

results_basic = run_tempo_basic(basic_settings)
settings = basic_settings

tim_file_path = joinpath(settings.work_dir, settings.tim_file)
backends = settings.backends

# Чтение данных из файла .tim
tim_file_data = readdlm(tim_file_path, String)
uncertainties = parse.(Float64, tim_file_data[3:end, 4])
residuals = readdlm(joinpath(settings.work_dir, "postfit.res"), Float64)[:, 2] .* 1e6


EFAC_max = 3.0
EQUAD_max = 200.0
lambda_max = 1000.0

EFAC_arr = LinRange(0, EFAC_max, 501)
EQUAD_arr = LinRange(0, EQUAD_max, 501)
lambda_arr = LinRange(0, lambda_max, 501)

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    adaptive_ks_fun = AdaptiveSparseGrid([0.0, 0.0], [EFAC_max, lambda_max], max_depth = 10, tol=1e-4) do (EFAC, lambda)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)
        ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
        return (ks_stat = KS_objective(res, unc, EFAC, lambda), 
        )
    end

    ks_arr = [EFAC == 0 ? Inf : adaptive_ks_fun([EFAC, lambda]) for EFAC in EFAC_arr, lambda in lambda_arr]

    ks_min_lambda = [minimum(ks_arr[:,i]) for i in eachindex(lambda_arr)]
    plot(lambda_arr, ks_min_lambda, label = backend)

end

legend()


fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function AD_objective(res, unc, EFAC, lambda)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)
        ad_test = OneSampleADTest(res_norm, Normal(0, 1))
        return ad_test.A² # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end


    # adaptive_ks_fun = AdaptiveSparseGrid([0.0, 0.0], [EFAC_max, lambda_max], max_depth = 10, tol=1e-4) do (EFAC, lambda)
    #     res_norm = res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)
    #     # ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
    #     # return (ks_stat = ks_test.δ, 
    #     # )
    #     ad_test = OneSampleADTest(res_norm, Normal(0, 1))
    #     println(EFAC, " ", lambda, " ", ad_test.A²)
    #     return (ks_stat = ad_test.A², 
    #     )
    # end

    # ks_arr = [EFAC == 0 ? Inf : adaptive_ks_fun([EFAC, lambda]) for EFAC in EFAC_arr, lambda in lambda_arr]

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : AD_objective(res, unc, EFAC, lambda) for EFAC in EFAC_arr, lambda in lambda_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_lambda = [minimum(ks_arr[:,i]) for i in eachindex(lambda_arr)]
    EQUAD_arr = [EFAC_arr[findmin(ks_arr[:,i])[2]] * lambda_arr[i] for i in eachindex(lambda_arr)]
    plot(EQUAD_arr, ks_min_lambda, ".", label = backend)

end

legend()




#-------------------------------------------------------------------------------------
# AD equad

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function AD_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ad_test = OneSampleADTest(res_norm, Normal(0, 1))
        return ad_test.A² # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end


    # adaptive_ks_fun = AdaptiveSparseGrid([0.0, 0.0], [EFAC_max, lambda_max], max_depth = 10, tol=1e-4) do (EFAC, lambda)
    #     res_norm = res ./ transform_uncertainty.(unc, EFAC, EFAC * lambda)
    #     # ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
    #     # return (ks_stat = ks_test.δ, 
    #     # )
    #     ad_test = OneSampleADTest(res_norm, Normal(0, 1))
    #     println(EFAC, " ", lambda, " ", ad_test.A²)
    #     return (ks_stat = ad_test.A², 
    #     )
    # end

    # ks_arr = [EFAC == 0 ? Inf : adaptive_ks_fun([EFAC, lambda]) for EFAC in EFAC_arr, lambda in lambda_arr]

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : AD_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EQUAD = [minimum(ks_arr[:,i]) for i in eachindex(EQUAD_arr)]
    plot(log10.(EQUAD_arr) .- 6.0, ks_min_EQUAD, ".", label = backend)

end
yscale("log")
xlabel("log10(EQUAD) - 6")
ylabel("Anderson Darling statistics")
legend()

#-------------------------------------------------------------------------------------
# KS equad

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function KS_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
        return sqrt(ks_test.n)*ks_test.δ # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : KS_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EQUAD = [minimum(ks_arr[:,i]) for i in eachindex(EQUAD_arr)]
    plot(log10.(EQUAD_arr) .- 6.0, ks_min_EQUAD, ".", label = backend)

end
yscale("log")
xlabel("log10(EQUAD) - 6")
ylabel("Kolmogorov Smirnov statistics")
legend()

#-------------------------------------------------------------------------------------
# JB equad

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function JB_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        jb_test = JarqueBeraTest(res_norm)
        return jb_test.JB # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : JB_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EQUAD = [minimum(ks_arr[:,i]) for i in eachindex(EQUAD_arr)]
    plot(log10.(EQUAD_arr) .- 6.0, ks_min_EQUAD, ".", label = backend)

end
yscale("log")
xlabel("log10(EQUAD) - 6")
ylabel("Jarque–Bera statistics")
legend()

#-------------------------------------------------------------------------------------
#  equad

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function loglikelihood_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        return loglikelihood(Normal(0,1), res_norm) # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? -Inf : loglikelihood_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EQUAD = [maximum(ks_arr[:,i]) for i in eachindex(EQUAD_arr)]
    plot(log10.(EQUAD_arr) .- 6.0, ks_min_EQUAD, ".", label = backend)

end
yscale("log")
xlabel("log10(EQUAD) - 6")
ylabel("loglikelihood ")
legend()

#-------------------------------------------------------------------------------------
# AD efac

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function AD_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ad_test = OneSampleADTest(res_norm, Normal(0, 1))
        return ad_test.A² # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : AD_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]
    plot(EFAC_arr, ks_min_EFAC, ".", label = backend)

end
yscale("log")
xlabel("EFAC")
ylabel("Anderson Darling statistics")
legend()


#-------------------------------------------------------------------------------------
# AD efac

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function AD_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ad_test = OneSampleADTest(res_norm, Normal(0, 1))
        return ad_test.A² # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : AD_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]
    plot(EFAC_arr, ks_min_EFAC, ".", label = backend)

end
yscale("log")
xlabel("EFAC")
ylabel("Anderson Darling statistics")
legend()

#-------------------------------------------------------------------------------------
# KS efac

fig, ax = subplots()

for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function KS_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
        return sqrt(ks_test.n)*ks_test.δ # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [EFAC == 0 ? Inf : KS_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]
    plot(EFAC_arr, ks_min_EFAC, ".", label = backend)

end
yscale("log")
xlabel("EFAC")
ylabel("KS statistics")
legend()


#-------------------------------------------------------------------------------------
# AD maps


for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function AD_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ad_test = OneSampleADTest(res_norm, Normal(0, 1))
        return ad_test.A² # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [AD_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]

    EFAC_min = EFAC_arr[findmin(ks_arr)[2][1]]
    EQUAD_min = EQUAD_arr[findmin(ks_arr)[2][2]]

    fig, ax = subplots()
    pcolormesh(EQUAD_arr, EFAC_arr, ks_arr, cmap="Blues_r", norm = matplotlib.colors.LogNorm(vmax=1e3))
    colorbar(label = "Anderson-Darling statistics")
    title(backend)
    xlabel("EQUAD")
    ylabel("EFAC")
    plot(EQUAD_min, EFAC_min, "x", color="red")

    legend()

end

#-------------------------------------------------------------------------------------
# KS maps


for backend in backends
    indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
    res = residuals[indices]
    unc = uncertainties[indices]

    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)

    function KS_objective(res, unc, EFAC, EQUAD)
        res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
        ks_test = ExactOneSampleKSTest(res_norm, Normal(0, 1))
        return sqrt(ks_test.n)*ks_test.δ # Минимизация p-значения или максимизация, в зависимости от вашей задачи
    end

    println(backend)

    ks_arr = [KS_objective(res, unc, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

    delete_nans(x) = isnan(x) ? Inf : x

    ks_arr .= delete_nans.(ks_arr)

    ks_min_EFAC = [minimum(ks_arr[i,:]) for i in eachindex(EFAC_arr)]

    EFAC_min = EFAC_arr[findmin(ks_arr)[2][1]]
    EQUAD_min = EQUAD_arr[findmin(ks_arr)[2][2]]

    fig, ax = subplots()
    pcolormesh(EQUAD_arr, EFAC_arr, ks_arr, cmap="Blues_r", norm = matplotlib.colors.LogNorm())
    colorbar(label = "KS statistics")
    title(backend)
    xlabel("EQUAD")
    ylabel("EFAC")
    plot(EQUAD_min, EFAC_min, "x", color="red")

    legend()

end
