#--------------------------------------------------------------------------------------------------------------
# Настройки для массива значений параметра (если у вас есть такой функционал)
struct ParameterSweepSettings
    parameter_name::String
    values::Vector{Float64}
    prior_mean::Union{Float64, Nothing}
    prior_std::Union{Float64, Nothing}
    prior_pdf::Union{Function, Nothing}
end

# Общий конструктор
ParameterSweepSettings(;parameter_name::String, values::Vector{Float64}=Float64[], prior_mean::Union{Float64, Nothing}=nothing, prior_std::Union{Float64, Nothing}=nothing, prior_pdf::Union{Function, Nothing}=nothing) = ParameterSweepSettings(parameter_name, values, prior_mean, prior_std, prior_pdf)

# Методы конструкторов для разных случаев
function ParameterSweepSettings(;parameter_name::String, prior_mean::Float64, prior_std::Float64)
    values = generate_values_from_prior(prior_mean, prior_std) # Генерация значений из нормального распределения
    return ParameterSweepSettings(parameter_name=parameter_name, values=values, prior_mean=prior_mean, prior_std=prior_std)
end

function ParameterSweepSettings(;parameter_name::String, prior_pdf::Function)
    values = generate_values_from_pdf(prior_pdf) # Генерация значений из PDF
    return ParameterSweepSettings(parameter_name=parameter_name, values=values, prior_pdf=prior_pdf)
end

# Пример использования
# settings1 = ParameterSweepSettings(parameter_name="XPBDOT", values=[1.0, 2.0, 3.0])
# settings2 = ParameterSweepSettings(parameter_name="XPBDOT", prior_mean=0.0, prior_std=1.0)
# settings3 = ParameterSweepSettings(parameter_name="XPBDOT", prior_pdf=my_pdf_function)

function run_tempo_parameter_sweep(
    basic_settings::BasicTempoSettings,
    global_iter_settings::GlobalIterationsSettings,
    parameter_sweep_settings::ParameterSweepSettings)
    # Результаты для каждого значения XPBDOT
    sweep_results = Dict()

    # Перебираем значения XPBDOT
    for value in parameter_sweep_settings.values
        # Создаем новый параметр для свипа
        sweep_parameter = GeneralTempoParameter(parameter_sweep_settings.parameter_name, value, flag=0)

        # Обновляем или добавляем параметр в настройки
        updated_tparams = update_or_add_tparam(copy(basic_settings.tparams), sweep_parameter)

        # Создаем новые базовые настройки с обновленным значением XPBDOT
        updated_basic_settings = BasicTempoSettings(
            work_dir = basic_settings.work_dir,
            version = basic_settings.version,
            par_file_init = basic_settings.par_file_init,
            tim_file = basic_settings.tim_file,
            flags = basic_settings.flags,
            tparams = updated_tparams,
            keys = basic_settings.keys
        )

        # Запускаем Tempo с обновленными базовыми настройками и текущими глобальными итерациями
        results_global_iters = run_tempo_global_iters(updated_basic_settings, global_iter_settings)

        # Сохраняем результаты
        sweep_results[value] = parsed_output
    end

    return sweep_results
end

