function calculate_EFACs_EQUADs(settings::BasicTempoSettings; time_start=-Inf, time_finish=+Inf)
    # Извлекаем пути к файлам и бэкенды из настроек
    tim_file_path = joinpath(settings.work_dir, settings.tim_file)
    backends = settings.backends

    # Чтение данных из файла .tim
    tim_file_data = readdlm(tim_file_path, String)
    times = parse.(Float64, tim_file_data[3:end, 3])
    uncertainties_orig = parse.(Float64, tim_file_data[3:end, 4])

    resdata = readdlm(joinpath(settings.work_dir, "residuals.dat"), Float64)

    residuals     = resdata[:,4] .* 1e6
    uncertainties = resdata[:,5] .* 1e6

    time_mask = (times .>= time_start) .& (times .<= time_finish)

    # Словари для результатов
    EFACs = Dict{String, Float64}()
    EQUADs = Dict{String, Float64}()
    log10EQUADs = Dict{String, Float64}()

    uncertainties_transformed = similar(residuals)
    residuals_shifted_norm    = similar(residuals)
    residuals_centered_norm   = similar(residuals)

    # Вычисление EFAC и EQUAD для каждого бэкенда
    for backend in backends
        
        indices_time = [i for i in 1:size(tim_file_data, 1) - 2 if (backend in tim_file_data[i+2, :] && time_mask[i] == 1)]  
        indices      = [i for i in 1:size(tim_file_data, 1) - 2 if (backend in tim_file_data[i+2, :])]  

        try
            EFAC, EQUAD, offset, objective = estimate_WhiteNoise_AD_with_offset(residuals[indices_time], uncertainties_orig[indices_time])
        
            uncertainties_transformed[indices] .= transform_uncertainties(uncertainties_orig[indices], EFAC, EQUAD)

            residuals_shifted_norm[indices]  .= (residuals[indices] .- offset) ./ uncertainties_transformed[indices]
            residuals_centered_norm[indices] .= (residuals[indices] .- weighted_mean(residuals[indices_time], uncertainties_transformed[indices_time])[1]) ./ uncertainties_transformed[indices]

            EFACs[backend] = EFAC
            EQUADs[backend] = EQUAD 
            if EQUAD > 0 
                log10EQUADs[backend] = log10(EQUAD) - 6.0
            else
                log10EQUADs[backend] = -10.0
            end
        catch error
            println(error)
            max_objective = Inf
            EFACs[backend] = 1.0
            EQUADs[backend] = 1e-4
            log10EQUADs[backend] = -10.0
        end
    end

    AD_objective = AD_objective_function(residuals_shifted_norm)
    chisqr_full  = sum(residuals_centered_norm .^ 2)

    # chisqr_full  = sum(((residuals .- weighted_mean(residuals, uncertainties)[1]) ./ uncertainties) .^ 2)

    # println(EFACs)
    # println(EQUADs)

    return EFACs, EQUADs, log10EQUADs, AD_objective, chisqr_full
end

function AD_objective_function_with_offset(residuals, uncertainties_transformed, offset)
	
	residuals_shifted = residuals .- offset

    residuals_shifted_norm = residuals_shifted ./ uncertainties_transformed

    ad_test = OneSampleADTest(residuals_shifted_norm, Normal(0, 1))

    return ad_test.A²
end

function estimate_WhiteNoise_AD_with_offset_old(residuals, uncertainties_orig; print_results = false, plot_results = false, backend = "")
	N = length(residuals)
	SD_mean = 1.0 / sqrt(N)
	SD_median = 1.0 / sqrt(2.0 * pi * N)
	SD_std = 1.0 / sqrt(2 * (N - 1))
	SD_skewness = sqrt(6.0 / N)
	SD_kurtosis = sqrt(24.0 / N)

	function EFAC_EQUAD_offset_AD_objective(args; residuals = residuals, uncertainties_orig = uncertainties_orig)
    	EFAC, EQUAD, offset = args

    	uncertainties_transformed = transform_uncertainty.(uncertainties_orig, EFAC, EQUAD)

    	AD_objective = AD_objective_function_with_offset(residuals, uncertainties_transformed, offset)

	   	return AD_objective
    end

    AD_result = optimize(EFAC_EQUAD_offset_AD_objective, [1.0, 0.0, 0.0])

    EFAC, EQUAD, offset = Optim.minimizer(AD_result)
    EFAC  = abs(EFAC)
    EQUAD = abs(EQUAD)


    function EFAC_EQUAD_offset_AD_std_objective(args; residuals = residuals, uncertainties_orig = uncertainties_orig)
    	EFAC, EQUAD, offset = args

    	uncertainties_transformed = transform_uncertainty.(uncertainties_orig, EFAC, EQUAD)

    	AD_objective = AD_objective_function_with_offset(residuals, uncertainties_transformed, offset)

    	residuals_norm = residuals ./ uncertainties_transformed
    	std_residuals_norm = std(residuals_norm)

    	std_objective = std_residuals_norm >= 1.0 ? 0.0 : (1.0 - std_residuals_norm) / 1e-6
    	return AD_objective + std_objective
    end

    AD_std_result = optimize(EFAC_EQUAD_offset_AD_std_objective, [EFAC, EQUAD, offset])

    EFAC, EQUAD, offset = Optim.minimizer(AD_std_result)
    EFAC  = abs(EFAC)
    EQUAD = abs(EQUAD)

    AD_std_objective = AD_std_result.minimum

    uncertainties_transformed = transform_uncertainty.(uncertainties_orig, EFAC, EQUAD)
    residuals_norm = residuals ./ uncertainties_transformed

    AD_objective = AD_objective_function_with_offset(residuals, uncertainties_transformed, offset)

    if print_results
    	println("backend: $backend")
    	println("	N_TOAs: $N")
		println("residuals stats:")
		println("	weighted mean:     $(weighted_mean(residuals, uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(residuals, uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(residuals, uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(residuals, uncertainties_transformed))")

		println("original uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_orig))")
		println("	weighted std:      $(weighted_std(uncertainties_orig))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_orig))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_orig))")

		println("transformed uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_transformed))")

		println("normalized residuals stats:")
		println("	mean:     $(weighted_mean(residuals_norm))")
		println("	std:      $(weighted_std(residuals_norm))")
		println("	skewness: $(weighted_skewness(residuals_norm))")
		println("	kurtosis: $(weighted_kurtosis(residuals_norm))")

		println("chi squared stats:")
		chi2, chi2r = chisq_stats(residuals, uncertainties_transformed)
		println("	chi2: $chi2")
		println("	chi2r: $chi2r")

		println("AD_objective:  $AD_objective")
		println("AD_std_objective:  $AD_std_objective")
		println("EFAC:   $EFAC")
		println("EQUAD:  $EQUAD")
		println("offset: $offset")

		println()




    	# println("Mean:     $(mean(residuals_norm)), $(mean(residuals_norm) / SD_mean)")
    	# println("Median:   $(median(residuals_norm)), $(median(residuals_norm) / SD_median)")
    	# println("Std:      $(std(residuals_norm)), $((std(residuals_norm) - 1.0) / SD_std)")
    	# println("Skewness: $(skewness(residuals_norm)), $(skewness(residuals_norm) / SD_std)")
    	# println("Kurtosis: $(kurtosis(residuals_norm)), $(kurtosis(residuals_norm) / SD_std)")
    	# println("EFAC = $EFAC, EQUAD = $EQUAD, offset = $offset, AD_objective = $(AD_objective)\n")
    end

    if plot_results
    	
		EFAC_arr  = collect(LinRange(0.01, ceil(std(residuals ./ uncertainties_orig)), 128))
		EQUAD_arr = collect(LinRange(0.0, ceil(std(residuals)), 128))
		AD_objective_arr = [AD_objective_function_offset(residuals, uncertainties_orig, EFAC, EQUAD) for EFAC in EFAC_arr, EQUAD in EQUAD_arr]
		std_arr = [weighted_std(residuals ./ transform_uncertainty.(uncertainties_orig, EFAC, EQUAD))[1] for EFAC in EFAC_arr, EQUAD in EQUAD_arr]

		fig, ax = subplots()
		imsh = ax.imshow(log10.(AD_objective_arr), extent=(EQUAD_arr[1], EQUAD_arr[end], EFAC_arr[1], EFAC_arr[end]), cmap="Blues_r", origin="lower", aspect="auto")
		cbar = colorbar(imsh)
		cbar.set_label(L"$\log_{10}\mathrm{AD_objective}$", fontsize=14)

		cs2 = ax.contour(EQUAD_arr, EFAC_arr, AD_objective_arr, levels = [0.5, 0.75, 1.0, 1.25], linestyles=["-", "--", "-.", ":"], colors="green")
		plot([], [], label="AD statistics", "-",  color="green")

		cs3 = ax.contour(EQUAD_arr, EFAC_arr, std_arr, levels = [1.0], linestyles=["-"], colors="violet")
		plot([], [], label="std", "-",  color="violet")

		plot(EQUAD, EFAC, "x", color="red")
		xlabel("EQUAD")
		ylabel("EFAC")
		title("$backend")
		legend()
		tight_layout()
	end

    return EFAC, EQUAD, AD_objective
end



transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)
kurtosis_function(res, unc, EFAC, EQUAD) = (EFAC == 0 && EQUAD == 0) ? Inf : kurtosis(res ./ transform_uncertainty.(unc, EFAC, EQUAD))

function AD_objective_function(residuals, uncertainties, EFAC, EQUAD; centered=true)
    residuals_norm = residuals ./ transform_uncertainty.(uncertainties, EFAC, EQUAD)
    residuals_norm_centered = residuals_norm .- centered .* mean(residuals_norm)
    ad_test = OneSampleADTest(residuals_norm_centered, Normal(0, 1))
    residuals_norm_std = std(residuals_norm_centered)
    return ad_test.A²
end

function estimate_WhiteNoise_AD(residuals, uncertainties; print_results = false)
	N = length(residuals)
	SD_mean = 1.0 / sqrt(N)
	SD_median = 1.0 / sqrt(2.0 * pi * N)
	SD_std = 1.0 / sqrt(2 * (N - 1))
	SD_skewness = sqrt(6.0 / N)
	SD_kurtosis = sqrt(24.0 / N)

	uncertainties_mean = mean(uncertainties)

    lambda_init = 0.1

    EQUAD_from_lambda(lambda, EFAC) = lambda * EFAC * uncertainties_mean

    function find_EFAC_normalizing_std(lambda; residuals = residuals, uncertainties = uncertainties)
		EFAC = std(residuals ./ transform_uncertainty.(uncertainties, 1.0, EQUAD_from_lambda(lambda, 1.0)))
   		return EFAC
	end

    function lambda_AD_objective(args; residuals = residuals, uncertainties = uncertainties)
    	lambda = args[1]
    	EFAC = find_EFAC_normalizing_std(lambda, residuals = residuals, uncertainties = uncertainties)
    	EQUAD = EQUAD_from_lambda(lambda, EFAC)
    	AD_objective = AD_objective_function(residuals, uncertainties, EFAC, EQUAD)
    	return AD_objective
    end

    AD_result = optimize(lambda_AD_objective, [lambda_init])

    lambda = abs(Optim.minimizer(AD_result)[1])
    EFAC = find_EFAC_normalizing_std(lambda)
    EQUAD = EQUAD_from_lambda(lambda, EFAC)
    AD_objective = AD_objective_function(residuals, uncertainties, EFAC, EQUAD)

    if print_results
    	uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
    	residuals_norm = residuals ./ uncertainties_transformed

    	println("Mean:     $(mean(residuals_norm)), $(mean(residuals_norm) / SD_mean)")
    	println("Median:   $(median(residuals_norm)), $(median(residuals_norm) / SD_median)")
    	println("Std:      $(std(residuals_norm)), $((std(residuals_norm) - 1.0) / SD_std)")
    	println("Skewness: $(skewness(residuals_norm)), $(skewness(residuals_norm) / SD_std)")
    	println("Kurtosis: $(kurtosis(residuals_norm)), $(kurtosis(residuals_norm) / SD_std)")
    	println("EFAC = $EFAC, EQUAD = $EQUAD, lambda = $lambda, AD_objective = $(AD_objective)\n")
    end

    return EFAC, EQUAD, lambda, AD_objective
end

function estimate_WhiteNoise_AD_outliers(residuals, uncertainties; max_iterations=30, AD_threshold=1.26, Z_threshold=3.0, print_results=false)
    current_residuals = copy(residuals)
    current_uncertainties = copy(uncertainties)
    AD_full_data = Inf

    if print_results
        println("Статистика неопределенностей: mean = $(mean(current_uncertainties)), std = $(std(current_uncertainties))")
    end

    for iter in 1:max_iterations
        # Оптимизируем EFAC и EQUAD на текущих данных
        EFAC, EQUAD, lambda, AD_value = estimate_WhiteNoise_AD(current_residuals, current_uncertainties, print_results=print_results)
        if iter == 1
            AD_full_data = AD_value
        end

        # Если AD-статистика уже <= AD_threshold, выходим
        if AD_value <= AD_threshold
            if print_results 
                println("AD уменьшен до $AD_value, завершение итераций.")
            end
            break
        end

        # Вычисляем текущее среднее и стандартное отклонение нормализованных резидуалов
        current_uncertainties_transformed = transform_uncertainty.(current_uncertainties, EFAC, EQUAD)
        current_residuals_norm = current_residuals ./ current_uncertainties_transformed
        residuals_mean = mean(current_residuals_norm)
        residuals_std = std(current_residuals_norm)

        # Вычисляем Z-score всех точек
        Z_scores = abs.(current_residuals_norm .- residuals_mean) ./ residuals_std

        # Находим индексы точек, у которых Z-score > Z_threshold
        candidate_indices = findall(Z_scores .> Z_threshold)

        if isempty(candidate_indices)
            if print_results
                println("Нет точек с Z-score > $Z_threshold. Завершение итераций.")
            end
            break
        end

        # Проверяем удаление только выбросов
        best_index = -1
        AD_best = AD_value
        best_Zscore = 0.0  # Запоминаем Z-score выброшенной точки

        for i in candidate_indices
            temp_residuals = deleteat!(copy(current_residuals), i)
            temp_uncertainties = deleteat!(copy(current_uncertainties), i)

            # Оптимизируем EFAC и EQUAD без этой точки
            EFAC_new, EQUAD_new, lambda_new, AD_new = estimate_WhiteNoise_AD(temp_residuals, temp_uncertainties)

            # Если удаление этой точки уменьшает AD, запоминаем её
            if AD_new < AD_best
                AD_best = AD_new
                best_index = i
                best_Zscore = Z_scores[i]  # Запоминаем Z-score выброшенной точки
            end
        end

        # Если ни одна точка не улучшает AD, выходим
        if best_index == -1
            if print_results
                println("Не найдено точек, удаление которых улучшает AD.")
            end
            break
        end

        # Удаляем лучшую точку и выводим её отклонение
        if print_results
            println("Удаление точки $best_index снижает AD с $AD_value до $AD_best (нормированное отклонение: $best_Zscore σ)")
        end
        deleteat!(current_residuals, best_index)
        deleteat!(current_uncertainties, best_index)
    end

    # Финальный пересчёт EFAC, EQUAD, lambda и AD-статистики на очищенных данных
    EFAC_final, EQUAD_final, lambda_final, AD_final = estimate_WhiteNoise_AD(current_residuals, current_uncertainties, print_results=print_results)

    # Выводим итоговые параметры
    if print_results
        println("\nФинальные параметры после удаления выбросов:")
        println("EFAC = $EFAC_final, EQUAD = $EQUAD_final, lambda = $lambda_final, AD_objective = $AD_final")
    end

    if print_results
        N = length(residuals)
        SD_mean = 1.0 / sqrt(N)
        SD_median = 1.0 / sqrt(2.0 * pi * N)
        SD_std = 1.0 / sqrt(2 * (N - 1))
        SD_skewness = sqrt(6.0 / N)
        SD_kurtosis = sqrt(24.0 / N)

        AD_objective = AD_objective_function(residuals, uncertainties, EFAC_final, EQUAD_final)
        uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC_final, EQUAD_final)
        residuals_norm = residuals ./ uncertainties_transformed

        println("Results on full dataset")
        println("Mean:     $(mean(residuals_norm)), $(mean(residuals_norm) / SD_mean)")
        println("Median:   $(median(residuals_norm)), $(median(residuals_norm) / SD_median)")
        println("Std:      $(std(residuals_norm)), $((std(residuals_norm) - 1.0) / SD_std)")
        println("Skewness: $(skewness(residuals_norm)), $(skewness(residuals_norm) / SD_std)")
        println("Kurtosis: $(kurtosis(residuals_norm)), $(kurtosis(residuals_norm) / SD_std)")
        println("EFAC = $EFAC_final, EQUAD = $EQUAD_final, lambda = $lambda_final, AD_objective = $(AD_objective)")
        println("AD_objective_full_data = $AD_full_data\n")
    end

    return EFAC_final, EQUAD_final, lambda_final, AD_full_data
end


function KS_objective(residuals, uncertainties, EFAC, EQUAD)
    residuals_norm = residuals ./ transform_uncertainty.(uncertainties, EFAC, EQUAD)
    ks_test = ApproximateOneSampleKSTest(residuals_norm, Normal(0, 1))
    residuals_norm_std = std(residuals_norm)
    return ks_test.δ + abs(residuals_norm_std - 1.0) # Минимизация статистики KS (максимального отклонения CDF)
end

function estimate_WhiteNoise_KS(residuals, uncertainties)
	sigma_cut = Inf

	uncertainties_mean = mean(uncertainties)

    EFAC   = 1.0
    lambda = 0.1
    EQUAD  = lambda * EFAC * uncertainties_mean

    mask = residuals .> -Inf
    
    for i in 1:10
        uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
        residuals_norm = residuals ./ uncertainties_transformed

        z_cut_level = std(residuals_norm[mask]) * sigma_cut
        mask .= abs.(residuals_norm .- mean(residuals_norm[mask])) .< z_cut_level

    	objective(EFAC_lambda) = KS_objective(residuals[mask], uncertainties[mask], EFAC_lambda[1], EFAC_lambda[2] * EFAC_lambda[1] * uncertainties_mean)

    	KS_result = optimize(objective, [EFAC, lambda])

    	EFAC, lambda = abs.(Optim.minimizer(KS_result))
    	EQUAD = lambda * EFAC * uncertainties_mean

    	# println("i = $i: EFAC = $EFAC, EQUAD = $EQUAD, obj = $(KS_result.minimum)")
    end

    # uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
    # residuals_norm = residuals ./ uncertainties_transformed

    # println("total # of TOAs: $(length(mask))")
    # println("# of discarded TOAs: $(length(mask) - sum(mask))")
    # println("All:    mean = $(mean(residuals_norm)), std = $(std(residuals_norm)), skew = $(skewness(residuals_norm)), kurt = $(kurtosis(residuals_norm))")
    # println("Masked: mean = $(mean(residuals_norm[mask])), std = $(std(residuals_norm[mask])), skew = $(skewness(residuals_norm[mask])), kurt = $(kurtosis(residuals_norm[mask]))")
    # println("EFAC = $EFAC, EQUAD = $EQUAD\n")

    return EFAC, EQUAD, mask
end

function moments_objective_trunc(res, unc, EFAC, EQUAD, sigma_cut)
    # Параметры обрезания
    lower_bound = - sigma_cut
    upper_bound = + sigma_cut
    dist = Normal(0, 1)

    # Функции для плотности и CDF
    phi(x) = pdf(dist, x)
    Phi(x) = cdf(dist, x)

    # Математическое ожидание
    mu_trunc = (phi(lower_bound) - phi(upper_bound)) / (Phi(upper_bound) - Phi(lower_bound))

    # Дисперсия
    sigma_trunc_sq = 1 + (lower_bound * phi(lower_bound) - upper_bound * phi(upper_bound)) / (Phi(upper_bound) - Phi(lower_bound)) - mu_trunc^2
    std_trunc = sqrt(sigma_trunc_sq)

    # Четвёртый момент
    function fourth_moment(mu_trunc, sigma_trunc_sq)
        integrand(x) = ((x - mu_trunc)^4 * phi(x)) / (Phi(upper_bound) - Phi(lower_bound))
        quadgk(integrand, lower_bound, upper_bound)[1]
    end
    fourth_mom = fourth_moment(mu_trunc, sigma_trunc_sq)

    # Куртосис
    kurtosis_trunc = fourth_mom / sigma_trunc_sq^2 - 3.0

    res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
    obj = (std_trunc / std(res_norm) - 1.0)^2 + (kurtosis(res_norm) - kurtosis_trunc)^2
    # println("EFAC = $EFAC, EQUAD = $EQUAD, obj = $obj, std = $(std(res_norm)), kurt = $(kurtosis(res_norm))")
    return obj # Минимизация p-значения или максимизация, в зависимости от вашей задачи
end

function estimate_WhiteNoise_moments_trunc(residuals, uncertainties)
    uncertainties_mean = mean(uncertainties)

    sigma_cut = 5.0

    EFAC = 1.0
    lambda = 0.1
    EQUAD = lambda * EFAC * uncertainties_mean
    mask = residuals .> -Inf
    
    for i in 1:10
        uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
        residuals_norm = residuals ./ uncertainties_transformed

        z_cut_level = std(residuals_norm[mask]) * sigma_cut
        mask .= abs.(residuals_norm .- mean(residuals_norm[mask])) .< z_cut_level

        objective(EFAC_lambda) = moments_objective_trunc(residuals[mask], uncertainties[mask], EFAC_lambda[1], EFAC_lambda[2] * EFAC_lambda[1] * uncertainties_mean, sigma_cut)

        moments_result = optimize(objective, [EFAC, lambda])

        EFAC, lambda = abs.(Optim.minimizer(moments_result))
        EQUAD = lambda * EFAC * uncertainties_mean

        # println("EFAC = $EFAC, lambda = $lambda, EQUAD = $EQUAD, obj = $(objective([EFAC, lambda])), #discarded = $(length(mask) - sum(mask))")
    end

    # uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
    # residuals_norm = residuals ./ uncertainties_transformed

    # println("total # of TOAs: $(length(mask))")
    # println("# of discarded TOAs: $(length(mask) - sum(mask))")
    # println("All:    mean = $(mean(residuals_norm)), std = $(std(residuals_norm)), skew = $(skewness(residuals_norm)), kurt = $(kurtosis(residuals_norm))")
    # println("Masked: mean = $(mean(residuals_norm[mask])), std = $(std(residuals_norm[mask])), skew = $(skewness(residuals_norm[mask])), kurt = $(kurtosis(residuals_norm[mask]))")
    # println("EFAC = $EFAC, EQUAD = $EQUAD\n")

    return EFAC, EQUAD
end

function moments_objective(res, unc, EFAC, EQUAD)
    res_norm = res ./ transform_uncertainty.(unc, EFAC, EQUAD)
    obj = (std(res_norm) - 1.0)^2 + (kurtosis(res_norm) - 0.0)^2
    return obj # Минимизация p-значения или максимизация, в зависимости от вашей задачи
end

function estimate_WhiteNoise_moments(residuals, uncertainties)

    uncertainties_mean = mean(uncertainties)

    EFAC = 1.0
    lambda = 0.1

    objective(EFAC_lambda) = moments_objective(residuals, uncertainties, EFAC_lambda[1], EFAC_lambda[2] * EFAC_lambda[1] * uncertainties_mean)

    for i in 1:3

        moments_result = optimize(objective, [EFAC, lambda])

        EFAC, lambda = abs.(Optim.minimizer(moments_result))
    end

    EQUAD = lambda * EFAC * uncertainties_mean

    # uncertainties_transformed = transform_uncertainty.(uncertainties, EFAC, EQUAD)
    # residuals_norm = residuals ./ uncertainties_transformed
    # println("All:    mean = $(mean(residuals_norm)), std = $(std(residuals_norm)), skew = $(skewness(residuals_norm)), kurt = $(kurtosis(residuals_norm))")
    # println("EFAC = $EFAC, EQUAD = $EQUAD\n")

    return EFAC, EQUAD
end

function update_EFACs_EQUADs_in_par_file!(par_file::TempoParFile, EFACs::Dict{String, Float64}, log10EQUADs::Dict{String, Float64})
    for (backend, EFAC) in EFACs
        efac_param = GeneralTempoParameter("TNEF -be $backend", EFAC)
        extend_par_file!(par_file, efac_param)
    end

    for (backend, log10EQUAD) in log10EQUADs
        equad_param = GeneralTempoParameter("TNEQ -be $backend", log10EQUAD)
        extend_par_file!(par_file, equad_param)
    end
end

# EFAC = 2.0; EQUAD = 0.0
# plot(residuals[indices], residuals[indices] ./ uncertainties[indices] ./ std(residuals[indices] ./ uncertainties[indices]), ".", color="blue")
# plot(residuals[indices], residuals[indices] ./ transform_uncertainty.(uncertainties[indices], 0.0, 100.0) ./ std(residuals[indices] ./ transform_uncertainty.(uncertainties[indices], 0.0, 100.0)), ".", color="red")
# plot(residuals[indices], residuals[indices] ./ transform_uncertainty.(uncertainties[indices], 3.0, 100.0) ./ std(residuals[indices] ./ transform_uncertainty.(uncertainties[indices], 3.0, 0.0)), ".", color="green")

# Функция estimate_WN должна быть определена для вычисления EFAC и EQUAD

# transform_uncertainties(uncertainty, EFAC, EQUAD) = sqrt(EFAC^2 * uncertainty^2 + EQUAD^2)

# function calculate_kurtosis(residuals, uncertainties, EFAC, EQUAD)
#     transformed_unc = transform_uncertainties.(uncertainties, EFAC, EQUAD)
#     return kurtosis(residuals ./ transformed_unc)
# end

# function optimize_lambda(residuals, uncertainties, initial_lambda)
#     lambda_function = lambda -> calculate_kurtosis(residuals, uncertainties, 1.0, lambda)
#     optimized_result = optimize(lambda_function, [initial_lambda])
#     return Optim.minimizer(optimized_result)[1]
# end

# function find_EFAC(residuals, uncertainties, lambda)
#     EFAC_function = EFAC -> std(residuals ./ transform_uncertainties.(uncertainties, EFAC, EFAC * lambda)) - 1.0
#     return find_zero(EFAC_function, 1.0)
# end


function split_timfile(settings::BasicTempoSettings; 
        time_start::Float64=-Inf, 
        time_finish::Float64=+Inf, 
        split_by_backend::Bool=true, 
        filter_by_time::Bool=true)
    # Извлекаем пути к файлам и бэкенды из настроек
    tim_file_path = joinpath(settings.work_dir, settings.tim_file)
    backends = settings.backends

    # Чтение данных из файла .tim
    tim_file_data = readdlm(tim_file_path, String)

    # Преобразуем времена в числовой формат
    times = parse.(Float64, tim_file_data[3:end, 3])

    # Если фильтрация по времени включена, создаем маску для отбора строк по времени
    time_mask = filter_by_time ? (times .>= time_start) .& (times .<= time_finish) : trues(length(times))

    # Функция для записи данных в файл
    function save_tim_file(file_name, indices)
        open(file_name, "w") do file
            # Сначала записываем первые две строки с настройками
            for i in 1:2
                println(file, join(tim_file_data[i, :], " "))
            end
            # Записываем строки, соответствующие индексам
            for i in indices
                println(file, join(tim_file_data[i, :], " "))
            end
        end
    end

    # Если нужно разбивать по бэкэндам
    if split_by_backend
        for backend in backends
        # Находим индексы строк, соответствующих бэкэнду и попадающих в диапазон времени
        indices = [i for i in 3:size(tim_file_data, 1) if backend in tim_file_data[i, :] && time_mask[i-2]]

        # Если диапазон времени не является бесконечным, добавим его к имени файла
        time_suffix = (time_start > -Inf || time_finish < Inf) ? "_$(time_start)_to_$(time_finish)" : ""
        backend_tim_file = joinpath(settings.work_dir, "$backend$time_suffix.tim")

        # Сохраняем файл
        save_tim_file(backend_tim_file, indices)
        println("Сохранен файл для бэкэнда $backend по пути: $backend_tim_file")
        end
    else
        # Если не нужно разбивать по бэкэндам, просто фильтруем строки по времени
        indices = [i for i in 3:size(tim_file_data, 1) if time_mask[i-2]]

        # Формируем имя файла, если диапазон времени задан
        time_suffix = (time_start > -Inf || time_finish < Inf) ? "_$(time_start)_to_$(time_finish)" : ""
        output_tim_file = joinpath(settings.work_dir, "filtered_data$time_suffix.tim")

        # Сохраняем файл
        save_tim_file(output_tim_file, indices)
        println("Сохранен файл с фильтрацией по времени по пути: $output_tim_file")
    end
end

#------------------------------------------------------------------------------------------------------------------------------------

function weighted_mean(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
	weights = unc .^ -2
	sum_weights = sum(weights)
	weighted_mean_res = sum(weights .* res) / sum_weights
	weighted_mean_unc = sqrt(1 / sum_weights)
	return weighted_mean_res, weighted_mean_unc
end

function weighted_std(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
	weights = unc .^ -2
	sum_weights = sum(weights)
	weighted_mean_res = sum(weights .* res) / sum_weights
	weighted_std_res = sqrt(sum(weights .* (res .- weighted_mean_res).^2) / sum_weights)
	N = length(res)
	N_eff = (sum_weights^2) / sum(weights.^2)
	weighted_std_res_unc = weighted_std_res / sqrt(2 * (N_eff - 1))
	return weighted_std_res, weighted_std_res_unc
end

function weighted_skewness(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)
    N = length(res)

    # Взвешенное среднее
    mean_w = sum(weights .* res) / sum_weights

    # Взвешенное стандартное отклонение
    std_w = sqrt(sum(weights .* (res .- mean_w).^2) / sum_weights)

    # Взвешенный skewness
    skew_numer = sum(weights .* (res .- mean_w).^3)
    skew_denom = sum_weights * std_w^3
    skew_w = skew_numer / skew_denom

    # Эффективное число наблюдений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка skewness
    skew_err = sqrt(6.0 / N_eff)

    return skew_w, skew_err
end

function weighted_kurtosis(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)
    N = length(res)

    # Взвешенное среднее
    mean_w = sum(weights .* res) / sum_weights

    # Взвешенное стандартное отклонение
    std_w = sqrt(sum(weights .* (res .- mean_w).^2) / sum_weights)

    # Взвешенный kurtosis (не эксцесс)
    kurt_numer = sum(weights .* (res .- mean_w).^4)
    kurt_denom = sum_weights * std_w^4
    kurt_w = kurt_numer / kurt_denom

    # Эффективное число наблюдений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка kurtosis
    kurt_err = sqrt(24.0 / N_eff)

    # Возвращаем также эксцесс
    excess_kurt_w = kurt_w - 3.0

    return excess_kurt_w, kurt_err
end

function weighted_rms(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)

    # RMS от нуля
    rms = sqrt(sum(weights .* res.^2) / sum_weights)

    # Эффективное число измерений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка RMS
    rms_err = rms / sqrt(2 * N_eff)

    return rms, rms_err
end

function chisq_stats(res::Vector{Float64}, unc::Vector{Float64}, dof::Int = length(res))
    # Предполагается, что unc уже включает EFAC и EQUAD
    chi2 = sum((res ./ unc).^2)
    red_chi2 = chi2 / dof
    return chi2, red_chi2
end

function transform_uncertainties!(uncertainties_transformed::Vector{Float64}, uncertainties_orig::Vector{Float64}, EFAC::Float64, EQUAD::Float64)
    @. uncertainties_transformed = sqrt(EFAC^2 * uncertainties_orig^2 + EQUAD^2)
    return uncertainties_transformed
end

function transform_uncertainties(uncertainties_orig::Vector{Float64}, EFAC::Float64, EQUAD::Float64)
	uncertainties_transformed = similar(uncertainties_orig)
	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
	return uncertainties_transformed
end


function AD_objective_function!(
    residuals_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64}
	)
    @. residuals_norm = residuals / uncertainties_transformed
    ad_test = OneSampleADTest(residuals_norm, Normal(0, 1))
    return ad_test.A²
end

function AD_objective_function(residuals_norm::Vector{Float64})
    ad_test = OneSampleADTest(residuals_norm, Normal(0, 1))
    return ad_test.A²
end

function AD_objective_function_with_offset!(
    residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64},
    offset::Float64
	)
    @. residuals_shifted_norm = (residuals - offset) / uncertainties_transformed
    ad_test = OneSampleADTest(residuals_shifted_norm, Normal(0, 1))
    return ad_test.A²
end


function AD_objective_function_fit_offset!(
	residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64}
	)

    function AD_objective_local(args)
    	offset = args[1]
        return AD_objective_function_with_offset!(
            residuals_shifted_norm,
            residuals,
            uncertainties_transformed,
            offset
        )
    end

    optim_result = optimize(AD_objective_local, [0.0])
    offset = optim_result.minimizer[1]
    AD_objective_val = AD_objective_local(offset)
    return AD_objective_val, offset
end

function find_EFAC_for_fixed_EQUAD_and_std!(
	residuals_norm::Vector{Float64},
	uncertainties_transformed::Vector{Float64},
	residuals::Vector{Float64},
	uncertainties_orig::Vector{Float64},
	EQUAD_fixed::Float64,
	std_res_norm_fixed::Float64
	)

    function std_res_norm_objective(EFAC::Float64)
        transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD_fixed)
        @. residuals_norm = residuals / uncertainties_transformed
        return std_res_norm_fixed / std(residuals_norm) - 1
    end

    std_residuals = std(residuals)

    if EQUAD_fixed > std_residuals / std_res_norm_fixed
        return NaN
    end

    mean_uncertainties_orig_squared = dot(uncertainties_orig, uncertainties_orig) / length(uncertainties_orig)

    EFAC_init = sqrt(((std_residuals / std_res_norm_fixed)^2 - EQUAD_fixed^2) / mean_uncertainties_orig_squared)

    EFAC = NaN

    try 
    	EFAC = find_zero(std_res_norm_objective, EFAC_init)
    catch error 
    	println(residuals)
    	println(uncertainties_orig)
    	println(EQUAD_fixed)
    	println(std_res_norm_fixed)
    	readline()
    end
    return abs(EFAC)
end

function estimate_WhiteNoise_AD_with_offset(residuals, uncertainties_orig; print_results = false, plot_results = false, backend = "", EFAC_init_in = nothing, EQUAD_init_in = nothing, offset_init_in = nothing)
	N_TOAs = length(residuals)
	SD_mean = 1.0 / sqrt(N_TOAs)
	SD_median = 1.0 / sqrt(2.0 * pi * N_TOAs)
	SD_std = 1.0 / sqrt(2 * (N_TOAs - 1))
	SD_skewness = sqrt(6.0 / N_TOAs)
	SD_kurtosis = sqrt(24.0 / N_TOAs)

	EFAC_max  = std(residuals ./ uncertainties_orig)
	EQUAD_max = std(residuals)

	uncertainties_transformed = similar(uncertainties_orig)
	residuals_norm            = similar(residuals)
	residuals_shifted_norm    = similar(residuals)
	
	function make_objective_function(
			residuals_shifted_norm::Vector{Float64},
			uncertainties_transformed::Vector{Float64},
			residuals::Vector{Float64},
			uncertainties_orig::Vector{Float64},
			)

    	function objective_fun(args)
        	EFAC, EQUAD, offset = args

        	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)

        	AD_objective = AD_objective_function_with_offset!(
            	residuals_shifted_norm,
            	residuals,
            	uncertainties_transformed,
            	offset
        	)

        	return AD_objective
    	end

    	return objective_fun
	end

	EFAC_EQUAD_offset_AD_objective = make_objective_function(residuals_shifted_norm, uncertainties_transformed, residuals, uncertainties_orig)


    if isnothing(EFAC_init_in) || isnothing(EQUAD_init_in) || isnothing(offset_init_in)
    	N_EQUADs = 11
    	N_std_res_norm = 11

    	EQUAD_init_arr = collect(LinRange(0.0, EQUAD_max, N_EQUADs))
    	std_res_norm_init_arr   = collect(LinRange(0.95, 1.05, N_std_res_norm))

    	EFAC_init_best, EQUAD_init_best, offset_init_best, AD_objective_init_best = NaN, NaN, NaN, Inf

    	for EQUAD_init in EQUAD_init_arr, std_res_norm_init in std_res_norm_init_arr 

    		EFAC_init = find_EFAC_for_fixed_EQUAD_and_std!(residuals_norm, uncertainties_transformed, residuals, uncertainties_orig, EQUAD_init, std_res_norm_init)

    		if isnan(EFAC_init)
    			continue
    		end

    		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC_init, EQUAD_init)

    		AD_objective_init, offset_init = AD_objective_function_fit_offset!(residuals_shifted_norm, residuals, uncertainties_transformed)

    		if AD_objective_init < AD_objective_init_best
    			AD_objective_init_best = AD_objective_init
    			EFAC_init_best         = EFAC_init
    			EQUAD_init_best        = EQUAD_init
    			offset_init_best       = offset_init
    		end
    	end

    	AD_result = optimize(EFAC_EQUAD_offset_AD_objective, [EFAC_init_best, EQUAD_init_best, offset_init_best])

    else
    	AD_result = optimize(EFAC_EQUAD_offset_AD_objective, [EFAC_init_in, EQUAD_init_in, offset_init_in])
    end

    EFAC_best, EQUAD_best, offset_best = Optim.minimizer(AD_result)
    EFAC_best  = abs(EFAC_best)
    EQUAD_best = abs(EQUAD_best)

    AD_objective_best = AD_result.minimum

    if print_results

    	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC_best, EQUAD_best)
    	@. residuals_norm = residuals / uncertainties_transformed

    	println("backend: $backend")
    	println("	N_TOAs: $N_TOAs")
		println("residuals stats:")
		println("	weighted mean:     $(weighted_mean(residuals, uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(residuals, uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(residuals, uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(residuals, uncertainties_transformed))")

		println("original uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_orig))")
		println("	weighted std:      $(weighted_std(uncertainties_orig))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_orig))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_orig))")

		println("transformed uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_transformed))")

		println("normalized residuals stats:")
		println("	mean:     $(weighted_mean(residuals_norm))")
		println("	std:      $(weighted_std(residuals_norm))")
		println("	skewness: $(weighted_skewness(residuals_norm))")
		println("	kurtosis: $(weighted_kurtosis(residuals_norm))")

		println("chi squared stats:")
		chi2, chi2r = chisq_stats(residuals, uncertainties_transformed)
		println("	chi2: $chi2")
		println("	chi2r: $chi2r")

		println("AD_objective_best:  $AD_objective_best")
		println("EFAC_best:   $EFAC_best")
		println("EQUAD_best:  $EQUAD_best")
		println("offset_best: $offset_best")

		println()

    end

    if plot_results
    	
		EFAC_arr  = collect(LinRange(0.01, ceil(EFAC_max), 128))
		EQUAD_arr = collect(LinRange(0.0, ceil(EQUAD_max), 128))
		AD_objective_arr = [
    		begin
        		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
        		AD_objective_function_fit_offset!(residuals_shifted_norm, residuals, uncertainties_transformed)[1]
    		end
    		for EFAC in EFAC_arr, EQUAD in EQUAD_arr
		]

		std_residuals_norm_arr = [
    		begin
        		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
        		@. residuals_norm = residuals / uncertainties_transformed
        		std(residuals_norm)
    		end
    		for EFAC in EFAC_arr, EQUAD in EQUAD_arr
		]

		fig, ax = subplots()
		imsh = ax.imshow(log10.(AD_objective_arr), extent=(EQUAD_arr[1], EQUAD_arr[end], EFAC_arr[1], EFAC_arr[end]), cmap="Blues_r", origin="lower", aspect="auto")
		cbar = colorbar(imsh)
		cbar.set_label(L"$\log_{10}\mathrm{AD_objective}$", fontsize=14)

		# cs2 = ax.contour(EQUAD_arr, EFAC_arr, AD_objective_arr, levels = [0.5, 0.75, 1.0, 1.25], linestyles=["-", "--", "-.", ":"], colors="green")
		# plot([], [], label="AD statistics", "-",  color="green")

		cs2 = ax.contour(EQUAD_arr, EFAC_arr, AD_objective_arr, levels = AD_objective_best .+ [0.1, 0.2, 0.3, 0.4], linestyles=["-", "--", "-.", ":"], colors="green")
		plot([], [], label="AD statistics", "-",  color="green")

		cs3 = ax.contour(EQUAD_arr, EFAC_arr, std_residuals_norm_arr, levels = [1.0], linestyles=["-"], colors="violet")
		plot([], [], label="std", "-",  color="violet")

		plot(EQUAD_best, EFAC_best, "x", color="red")

		AD_objective_grid, min_ind = findmin(AD_objective_arr)
		plot(EQUAD_arr[min_ind[2]], EFAC_arr[min_ind[1]], "x", color="black")

		println("AD_objective: optimized = $AD_objective_best, grid = $AD_objective_grid, delta = $(AD_objective_best - AD_objective_grid)")
		println("std: optimized = $(std(residuals_norm)), grid = $(std_residuals_norm_arr[min_ind])")

		xlabel("EQUAD")
		ylabel("EFAC")
		title("$backend")
		legend()
		tight_layout()

	end

    return EFAC_best, EQUAD_best, offset_best, AD_objective_best
end