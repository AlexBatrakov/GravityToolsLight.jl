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
    parameter_sweep_settings::ParameterSweepSettings)
    # Результаты для каждого значения XPBDOT
    sweep_results = GeneralTempoResult[]

    println("sweep parameter = $(parameter_sweep_settings.parameter_name)")

    # Перебираем значения XPBDOT
    for value in parameter_sweep_settings.values
        # Создаем новый параметр для свипа
        sweep_parameter = GeneralTempoParameter(parameter_sweep_settings.parameter_name, value, flag=0)

        # Обновляем или добавляем параметр в настройки
        updated_tparams = update_or_add_tparam(copy(basic_settings.tparams), sweep_parameter)
        print("value = $value, ")

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
        results_basic = run_tempo_basic(updated_basic_settings)
        chisqr = (results_basic.last_internal_iteration.result.chisqr)
        println("chisqr = $(chisqr), ") 

        # Сохраняем результаты
        push!(sweep_results, results_basic)
    end

    return sweep_results
end


function run_tempo_parameter_sweep(
    basic_settings::BasicTempoSettings,
    global_iters_settings::GlobalIterationsSettings,
    parameter_sweep_settings::ParameterSweepSettings)
    # Результаты для каждого значения XPBDOT
    sweep_results = GeneralTempoResult[]

    println("sweep parameter = $(parameter_sweep_settings.parameter_name)")

    # Перебираем значения XPBDOT
    for value in parameter_sweep_settings.values
        # Создаем новый параметр для свипа
        sweep_parameter = GeneralTempoParameter(parameter_sweep_settings.parameter_name, value, flag=0)

        # Обновляем или добавляем параметр в настройки
        updated_tparams = update_or_add_tparam(copy(basic_settings.tparams), sweep_parameter)
        println("value = $value ")

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
        results_global_iters = run_tempo_global_iters(updated_basic_settings, global_iters_settings)

        results = results_global_iters.last_internal_iteration.result

        println("chisqr = $(results.chisqr), chisqr/nfree = $(results.chisqr_nfree), TRES = $(results.TRES.post_fit), Pre/post = $(results.pre_post)")

        # Сохраняем результаты
        push!(sweep_results, results_global_iters)
    end

    return sweep_results
end
