function calculate_EFACs_EQUADs(settings::BasicTempoSettings)
    # Извлекаем пути к файлам и бэкенды из настроек
    tim_file_path = joinpath(settings.work_dir, settings.tim_file)
    backends = settings.backends

    # Чтение данных из файла .tim
    tim_file_data = readdlm(tim_file_path, String)
    uncertainties = parse.(Float64, tim_file_data[3:end, 4])
    residuals = readdlm(joinpath(settings.work_dir, "postfit.res"), Float64)[:, 2] .* 1e6

    # Словари для результатов
    EFACs = Dict{String, Float64}()
    EQUADs = Dict{String, Float64}()
    log10EQUADs = Dict{String, Float64}()

    # Вычисление EFAC и EQUAD для каждого бэкенда
    for backend in backends
        indices = [i for i in 1:size(tim_file_data, 1)-2 if backend in tim_file_data[i+2, :]]  
        
        try
        EFAC, EQUAD, indices_cut = GravityToolsLight.estimate_WhiteNoise(residuals, uncertainties, indices)
        
        EFACs[backend] = EFAC
        EQUADs[backend] = EQUAD
        if EQUAD > 0 
            log10EQUADs[backend] = log10(EQUAD) - 6.0
        end
        catch error
        
        end
    end

    return EFACs, EQUADs, log10EQUADs
end

function estimate_WhiteNoise(residuals, uncertainties, indices)
    # Определение функций для оптимизации
    transform_uncertainty(unc, EFAC, EQUAD) = sqrt(EFAC^2 * unc^2 + EQUAD^2)
    kurtosis_function(res, unc, EFAC, EQUAD) = (EFAC == 0 && EQUAD == 0) ? Inf : kurtosis(res ./ transform_uncertainty.(unc, EFAC, EQUAD))

    # Минимизация куртозиса для нахождения lambda
    lambda_function(lambda) = kurtosis_function(residuals[indices], uncertainties[indices], 1.0, lambda)
    lambda_opt = optimize(lambda_function, [0.0])
    lambda = max(Optim.minimizer(lambda_opt)[1], 0.0)

    # Нахождение EFAC
    EFAC_function(EFAC) = std(residuals[indices] ./ transform_uncertainty.(uncertainties[indices], EFAC, EFAC * lambda)) - 1.0
    EFAC = find_zero(EFAC_function, 1.0)
    EQUAD = lambda * EFAC


    unc_trans = transform_uncertainty.(uncertainties[indices], EFAC, EQUAD)

    hist(residuals[indices] ./ unc_trans)

    indices_cut = indices

    return EFAC, EQUAD, indices_cut
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